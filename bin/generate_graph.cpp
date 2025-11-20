#include <htslib/sam.h>
#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <getopt.h>
#include <cctype>
#include <iomanip>
#include <set>

// ===============================================================================
// 参数配置
// ===============================================================================

// 定义 contig 的端点区域：距离端点多远算作 START/END 区域
static int MAX_END = 1000;

// 质量过滤参数
static int MIN_MAPQ = 20;           // 最小比对质量
static int MAX_NM = 5;              // 最大错配数
static double MIN_MATCH_FRAC = 0.70; // 最小匹配比例（匹配碱基数/参考长度）

// 是否启用配对末端证据
static bool ENABLE_PAIRED = true;

// 最大跨度比例：read 距离端点的距离不能超过 contig 长度的这个比例
static double MAX_SPAN_FRAC = 0.80;

// 是否输出两种顺序的连接（A->B 和 B->A）
static bool OUTPUT_BOTH_ORDER = false;

// 测序文库类型：FR (forward-reverse), RF (reverse-forward), FF (forward-forward)
static std::string LIB_TYPE = "FR";

// 输出过滤参数
static int MIN_COUNT = 1;           // 最小支持 reads 数
static double MIN_SCORE = 0.0;      // 最小累积分数

// 调试模式开关
static bool DEBUG_MODE = false;

// ===============================================================================
// Contig 区域定义
// ===============================================================================

// 定义 contig 的三个区域：起始端、末尾端、中间
enum ContigRegion { START=0, END=1, MIDDLE=2 };

/**
 * 判断给定位置在 contig 的哪个区域
 */
static inline ContigRegion getContigRegion(int pos1, int contigLen){
    int pref = std::min(MAX_END, contigLen/2);
    int suff = std::max(contigLen - MAX_END, contigLen/2);
    if(pos1 <= pref) return START;
    if(pos1 > suff) return END;
    return MIDDLE;
}

/**
 * 计算位置到 contig 起始端的距离
 */
static inline int distToStart(int pos){
    return std::max(0, pos - 1);
}

/**
 * 计算位置到 contig 末尾端的距离
 */
static inline int distToEnd(int pos, int L){
    return std::max(0, L - pos);
}

/**
 * 翻转区域：START <-> END，MIDDLE 保持不变
 */
static inline ContigRegion flipRegion(ContigRegion r){
    if(r == START) return END;
    if(r == END) return START;
    return MIDDLE;
}

/**
 * 去除字符串首尾空白
 */
static inline void trim(std::string& s){
    size_t i = 0;
    while(i < s.size() && std::isspace((unsigned char)s[i])) ++i;
    size_t j = s.size();
    while(j > i && std::isspace((unsigned char)s[j-1])) --j;
    s = s.substr(i, j-i);
}

// ===============================================================================
// FastG 连接关系数据结构
// ===============================================================================

struct ContigPair {
    std::string ref1;
    std::string ref2;
    char orient1;  // '+' or '-'
    char orient2;  // '+' or '-'

    bool operator<(const ContigPair& other) const {
        if (ref1 != other.ref1) return ref1 < other.ref1;
        if (ref2 != other.ref2) return ref2 < other.ref2;
        if (orient1 != other.orient1) return orient1 < other.orient1;
        return orient2 < other.orient2;
    }
};

/**
 * 解析 FastG FAI 文件，获取预期的 contig 连接关系
 */
static std::set<ContigPair> parseFastgFile(const std::string& filename) {
    std::ifstream inFile(filename);
    std::string line;
    std::set<ContigPair> contigPairs;

    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        std::string fullName, contigName, linkedContig;
        std::getline(ss, fullName, ';');
        std::stringstream ss2(fullName);
        std::getline(ss2, contigName, ':');

        // 判断当前 contig 的方向
        bool contigReversed = false;
        if (!contigName.empty() && contigName.back() == '\'') {
            contigReversed = true;
            contigName.pop_back();
        }

        // 读取连接的 contigs
        while (std::getline(ss2, linkedContig, ',')) {
            if (linkedContig.empty()) continue;

            // 判断连接 contig 的方向
            bool linkedReversed = false;
            if (linkedContig.back() == '\'') {
                linkedReversed = true;
                linkedContig.pop_back();
            }

            // 计算方向
            char orient1, orient2;
            if (!contigReversed) {
                orient1 = '+';
                orient2 = linkedReversed ? '-' : '+';
            } else {
                orient1 = '-';
                orient2 = linkedReversed ? '+' : '-';
            }

            ContigPair pair{contigName, linkedContig, orient1, orient2};
            orient1 = (orient1 == '+') ? '-' : '+';
            orient2 = (orient2 == '+') ? '-' : '+';
            ContigPair con_pair{linkedContig, contigName, orient1, orient2};
            contigPairs.insert(pair);
            contigPairs.insert(con_pair);
        }
    }

    return contigPairs;
}

// ===============================================================================
// SA (Supplementary Alignment) 标签解析
// ===============================================================================

struct SAItem {
    std::string rname;
    int pos = -1;
    bool isRev = false;
    std::string cigar;
    int mapq = 0;
    int nm = 0;
    bool ok = false;
};

static SAItem parseSAItem(const std::string& item){
    SAItem it;
    std::stringstream ss(item);
    std::string rname, pos, strand, cigar, mapq, nm;

    if(!std::getline(ss, rname, ',') || !std::getline(ss, pos, ',') ||
       !std::getline(ss, strand, ',') || !std::getline(ss, cigar, ',') ||
       !std::getline(ss, mapq, ',') || !std::getline(ss, nm, ','))
        return it;

    trim(rname); trim(pos); trim(strand); trim(cigar); trim(mapq); trim(nm);
    if(rname.empty() || pos.empty()) return it;

    it.rname = rname;
    it.pos = std::atoi(pos.c_str());
    it.isRev = (strand == "-");
    it.cigar = cigar;
    it.mapq = std::atoi(mapq.c_str());
    it.nm = std::atoi(nm.c_str());
    it.ok = true;
    return it;
}

// ===============================================================================
// CIGAR 字符串处理
// ===============================================================================

static int cigarRefLen(const std::string& cigar){
    if(cigar.empty()) return 0;
    int n = 0, total = 0;
    for(char c: cigar){
        if(std::isdigit((unsigned char)c))
            n = n * 10 + (c - '0');
        else {
            if(c == 'M' || c == '=' || c == 'X' || c == 'D' || c == 'N')
                total += n;
            n = 0;
        }
    }
    return total;
}

static int cigarMatchLen(const std::string& cigar){
    if(cigar.empty()) return 0;
    int n = 0, total = 0;
    for(char c: cigar){
        if(std::isdigit((unsigned char)c))
            n = n * 10 + (c - '0');
        else {
            if(c == 'M' || c == '=' || c == 'X')
                total += n;
            n = 0;
        }
    }
    return total;
}

// ===============================================================================
// 过滤和评分函数
// ===============================================================================

static inline bool passMapqNm(int mapq, int nm){
    return mapq >= MIN_MAPQ && nm <= MAX_NM;
}

static inline bool passMatchFrac(int matchLen, int refLen){
    if(refLen <= 0) return false;
    return (double)matchLen / refLen >= MIN_MATCH_FRAC;
}

static inline double endWeight(int d1, int d2){
    double lambda = std::max(50.0, (double)MAX_END / 2.0);
    double w1 = std::exp(-(double)d1 / lambda);
    double w2 = std::exp(-(double)d2 / lambda);
    return w1 * w2;
}

// ===============================================================================
// 数据结构
// ===============================================================================

struct Evidence {
    std::string A, B;
    int LA = 0, LB = 0;
    int posA = 0, posB = 0;
    bool readOnAisRev = false;
    bool readOnBisRev = false;
    ContigRegion regA = MIDDLE;
    ContigRegion regB = MIDDLE;
    int mapqA = 0, nmA = 0;
    int mapqB = 0, nmB = 0;
    std::string readName;
    uint16_t readFlag = 0;
};

struct LayoutKey {
    std::string left, right;
    char oL, oR;

    LayoutKey(){}
    LayoutKey(const std::string& l, char ol, const std::string& r, char orr)
            : left(l), right(r), oL(ol), oR(orr) {}

    bool operator<(const LayoutKey& o) const {
        if(left != o.left) return left < o.left;
        if(right != o.right) return right < o.right;
        if(oL != o.oL) return oL < o.oL;
        return oR < o.oR;
    }
};

struct ReadInfo {
    std::string name;
    uint16_t flag;
    ReadInfo(const std::string& n, uint16_t f) : name(n), flag(f) {}
};

struct AggStats {
    int supplementCount = 0;        // supplement reads 计数（支持 fastg）
    int spanCount = 0;              // span read pairs 计数（支持 fastg）
    int supplementCountNoFastg = 0; // supplement reads 计数（不支持 fastg）
    int spanCountNoFastg = 0;       // span read pairs 计数（不支持 fastg）
    std::vector<ReadInfo> supportingReads;
};

// ===============================================================================
// 计算距离
// ===============================================================================

static void nearEndDistances(ContigRegion regL, int posL, int LL, char oL,
                             ContigRegion regR, int posR, int LR, char oR,
                             int& dNearL, int& dNearR){
    ContigRegion gRegL = (oL == '-') ? flipRegion(regL) : regL;
    ContigRegion gRegR = (oR == '-') ? flipRegion(regR) : regR;

    dNearL = (gRegL == START) ? distToStart(posL) : distToEnd(posL, LL);
    dNearR = (gRegR == START) ? distToStart(posR) : distToEnd(posR, LR);
}

// ===============================================================================
// CIGAR 区间解析
// ===============================================================================

struct ReadInterval {
    int start = 0;
    int end = 0;
    int len = 0;
    int softClipStart = 0;
    int softClipEnd = 0;
};

static ReadInterval parseCigarReadInterval(const std::string& cigar, bool isRev, int readLen = 0) {
    ReadInterval result;
    if(cigar.empty()) return result;

    std::vector<std::pair<int, char>> ops;
    int n = 0;
    for(char c : cigar) {
        if(std::isdigit((unsigned char)c)) {
            n = n * 10 + (c - '0');
        } else {
            if(n > 0) ops.push_back({n, c});
            n = 0;
        }
    }

    int readConsumed = 0;
    int softClipStart = 0;
    int softClipEnd = 0;

    if(!ops.empty() && ops[0].second == 'S') {
        softClipStart = ops[0].first;
    }

    if(ops.size() > 1 && ops.back().second == 'S') {
        softClipEnd = ops.back().first;
    }

    for(const auto& op : ops) {
        char c = op.second;
        int len = op.first;
        if(c == 'M' || c == 'I' || c == 'S' || c == '=' || c == 'X') {
            readConsumed += len;
        }
    }

    result.softClipStart = softClipStart;
    result.softClipEnd = softClipEnd;
    result.len = readConsumed;

    if(!isRev) {
        result.start = softClipStart + 1;
        result.end = readConsumed - softClipEnd;
    } else {
        if(readLen > 0) {
            result.start = readLen - (readConsumed - softClipEnd) + 1;
            result.end = readLen - softClipStart;
        } else {
            result.start = softClipStart + 1;
            result.end = readConsumed - softClipEnd;
        }
    }

    return result;
}

static int getReadLength(bam1_t* b) {
    int len = 0;
    uint32_t* cigar = bam_get_cigar(b);
    for(int i = 0; i < b->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int oplen = bam_cigar_oplen(cigar[i]);
        if(op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP ||
           op == BAM_CEQUAL || op == BAM_CDIFF) {
            len += oplen;
        }
    }
    return len;
}

// ===============================================================================
// 拼接检查
// ===============================================================================

static bool canStitchReadIntervals(const ReadInterval& interval1, const ReadInterval& interval2,
                                   int maxGap, int maxOverlap, bool& first1) {
    if(interval1.end <= interval2.start) {
        int gap = interval2.start - interval1.end - 1;
        if(gap <= maxGap) {
            first1 = true;
            return true;
        }
    }

    if(interval2.end <= interval1.start) {
        int gap = interval1.start - interval2.end - 1;
        if(gap <= maxGap) {
            first1 = false;
            return true;
        }
    }

    if(interval1.start <= interval2.end && interval2.start <= interval1.end) {
        int overlap = std::min(interval1.end, interval2.end) - std::max(interval1.start, interval2.start) + 1;
        if(overlap <= maxOverlap) {
            first1 = (interval1.start <= interval2.start);
            return true;
        }
    }

    return false;
}

// ===============================================================================
// 计算分数
// ===============================================================================

static bool computeLayoutScore(const Evidence& ev, bool leftIsA, char oL, char oR,
                               double& scoreOut) {
    int LL = leftIsA ? ev.LA : ev.LB;
    int LR = leftIsA ? ev.LB : ev.LA;
    int posL = leftIsA ? ev.posA : ev.posB;
    int posR = leftIsA ? ev.posB : ev.posA;
    ContigRegion regL = leftIsA ? ev.regA : ev.regB;
    ContigRegion regR = leftIsA ? ev.regB : ev.regA;
    int mapqL = leftIsA ? ev.mapqA : ev.mapqB;
    int nmL = leftIsA ? ev.nmA : ev.nmB;
    int mapqR = leftIsA ? ev.mapqB : ev.mapqA;
    int nmR = leftIsA ? ev.nmB : ev.nmA;

    int dNearL = 0, dNearR = 0;
    nearEndDistances(regL, posL, LL, oL, regR, posR, LR, oR, dNearL, dNearR);

    double w_end = endWeight(dNearL, dNearR);
    double w_qualL = std::min(1.0, (double)mapqL / 60.0) * (1.0 / (1.0 + 0.2 * std::max(0, nmL)));
    double w_qualR = std::min(1.0, (double)mapqR / 60.0) * (1.0 / (1.0 + 0.2 * std::max(0, nmR)));

    scoreOut = w_end * w_qualL * w_qualR;

    if(DEBUG_MODE) {
        std::cerr << "Score calculation: w_end=" << w_end
                  << " w_qualL=" << w_qualL << " w_qualR=" << w_qualR
                  << " total=" << scoreOut << "\n";
    }

    return (scoreOut > 0.0);
}

// ===============================================================================
// 配对末端检查
// ===============================================================================

static bool checkPairedEndLayout(int pos1, bool rev1, ContigRegion reg1, int L1,
                                 int pos2, bool rev2, ContigRegion reg2, int L2,
                                 char oL, char oR, bool first1) {
    bool revL, revR;
    ContigRegion regL, regR;
    int posL, posR, LL, LR;

    if(first1) {
        revL = rev1; revR = rev2;
        regL = reg1; regR = reg2;
        posL = pos1; posR = pos2;
        LL = L1; LR = L2;
    } else {
        revL = rev2; revR = rev1;
        regL = reg2; regR = reg1;
        posL = pos2; posR = pos1;
        LL = L2; LR = L1;
    }

    bool readIsForwardL = (oL == '-') ? revL : !revL;
    bool readIsForwardR = (oR == '-') ? revR : !revR;

    if(!readIsForwardL || readIsForwardR) return false;

    if(regL == MIDDLE || regR == MIDDLE) return false;

    ContigRegion leftContigPhysicalRightEnd = (oL == '+') ? END : START;
    ContigRegion rightContigPhysicalLeftEnd = (oR == '+') ? START : END;

    if(regL != leftContigPhysicalRightEnd) return false;
    if(regR != rightContigPhysicalLeftEnd) return false;

    int distL = (regL == START) ? distToStart(posL) : distToEnd(posL, LL);
    int distR = (regR == START) ? distToStart(posR) : distToEnd(posR, LR);

    double fracL = (LL > 0) ? ((double)distL / LL) : 1.0;
    double fracR = (LR > 0) ? ((double)distR / LR) : 1.0;

    if(fracL > MAX_SPAN_FRAC || fracR > MAX_SPAN_FRAC) return false;

    return true;
}

// ===============================================================================
// Split-read 检查
// ===============================================================================

static bool checkSplitReadLayout(int pos1, bool rev1, ContigRegion reg1, int L1,
                                 int pos2, bool rev2, ContigRegion reg2, int L2,
                                 char oL, char oR, bool first1) {
    bool revL_original, revR_original;
    ContigRegion regL_original, regR_original;

    if(first1) {
        revL_original = rev1; revR_original = rev2;
        regL_original = reg1; regR_original = reg2;
    } else {
        revL_original = rev2; revR_original = rev1;
        regL_original = reg2; regR_original = reg1;
    }

    bool readIsForwardL = (oL == '-') ? revL_original : !revL_original;
    bool readIsForwardR = (oR == '-') ? revR_original : !revR_original;

    if(!readIsForwardL || !readIsForwardR) return false;

    if(regL_original == MIDDLE || regR_original == MIDDLE) return false;

    ContigRegion leftContigPhysicalRightEnd = (oL == '+') ? END : START;
    ContigRegion rightContigPhysicalLeftEnd = (oR == '+') ? START : END;

    if(regL_original != leftContigPhysicalRightEnd) return false;
    if(regR_original != rightContigPhysicalLeftEnd) return false;

    return true;
}

// ===============================================================================
// 命令行帮助
// ===============================================================================

static void usage(const char* prog){
    std::cerr << "Usage: " << prog << " [options] <BAM> <FASTG_FAI> <OUT> <AverageDepth>\n";
    std::cerr << "Options:\n";
    std::cerr << "  -e <int>      MAX_END (default: 300)\n";
    std::cerr << "  -q <int>      MIN_MAPQ (default: 0)\n";
    std::cerr << "  -n <int>      MAX_NM (default: 10)\n";
    std::cerr << "  -p <double>   MIN_MATCH_FRAC (default: 0)\n";
    std::cerr << "  -P <0/1>      Enable paired-end evidence (default: 1)\n";
    std::cerr << "  --max-span-frac <double>  (default: 0.80)\n";
    std::cerr << "  --both-order <0/1>        Output both orders (default: 0)\n";
    std::cerr << "  --lib <FR|RF|FF>          Library type (default: FR)\n";
    std::cerr << "  --min-count <int>         Minimum supporting reads (default: 1)\n";
    std::cerr << "  --min-score <double>      Minimum accumulated score (default: 0.0)\n";
    std::cerr << "  --debug                   Enable debug output\n";
}

// ===============================================================================
// 主程序
// ===============================================================================

int main(int argc, char** argv){
    static struct option long_opts[] = {
            {"max-span-frac", required_argument, 0, 1000},
            {"both-order", required_argument, 0, 1001},
            {"lib", required_argument, 0, 1002},
            {"min-count", required_argument, 0, 1003},
            {"min-score", required_argument, 0, 1004},
            {"debug", no_argument, 0, 1005},
            {0,0,0,0}
    };
    int long_idx = 0;

    int opt;
    while((opt = getopt_long(argc, argv, "e:q:n:p:P:", long_opts, &long_idx)) != -1){
        switch(opt){
            case 'e': MAX_END = std::max(1, std::atoi(optarg)); break;
            case 'q': MIN_MAPQ = std::max(0, std::atoi(optarg)); break;
            case 'n': MAX_NM = std::max(0, std::atoi(optarg)); break;
            case 'p': MIN_MATCH_FRAC = std::max(0.0, std::min(1.0, std::atof(optarg))); break;
            case 'P': ENABLE_PAIRED = (std::atoi(optarg) != 0); break;
            case 1000: MAX_SPAN_FRAC = std::min(0.99, std::max(0.1, std::atof(optarg))); break;
            case 1001: OUTPUT_BOTH_ORDER = (std::atoi(optarg) != 0); break;
            case 1002: {
                std::string v = optarg;
                if(v=="FR"||v=="RF"||v=="FF") LIB_TYPE=v;
                else { std::cerr << "Unknown --lib " << v << ", fallback FR\n"; LIB_TYPE="FR"; }
                break;
            }
            case 1003: MIN_COUNT = std::max(1, std::atoi(optarg)); break;
            case 1004: MIN_SCORE = std::max(0.0, std::atof(optarg)); break;
            case 1005: DEBUG_MODE = true; break;
            default: usage(argv[0]); return 1;
        }
    }

    if(argc - optind < 4){
        usage(argv[0]);
        return 1;
    }

    const char* bamPath = argv[optind];
    const char* fastgFaiPath = argv[optind+1];
    const char* outPath = argv[optind+2];
    double avgDepth = std::atof(argv[optind+3]);

    // 解析 FastG FAI 文件
    std::set<ContigPair> fastgPairs = parseFastgFile(fastgFaiPath);

    if(DEBUG_MODE) {
        std::cerr << "Loaded " << fastgPairs.size() << " expected connections from FastG\n";
    }

    samFile* in = sam_open(bamPath, "r");
    if(!in){
        std::cerr << "Failed to open BAM " << bamPath << "\n";
        return 1;
    }

    bam_hdr_t* hdr = sam_hdr_read(in);
    if(!hdr){
        std::cerr << "Failed to read BAM header\n";
        sam_close(in);
        return 1;
    }

    std::unordered_map<std::string, int> name_to_tid;
    for(int i = 0; i < hdr->n_targets; i++){
        name_to_tid[hdr->target_name[i]] = i;
    }

    bam1_t* b = bam_init1();

    std::unordered_map<std::string, double> refConsumed;
    std::map<LayoutKey, AggStats> agg;

    std::unordered_set<std::string> processedSupplementReads;
    std::unordered_set<std::string> processedPairedReads;

    //auto consumeRefLen = [&](const char* rname, int n_cigar, uint32_t* cigar){
    //    int L = bam_cigar2rlen(n_cigar, cigar);
    //    if(L > 0) refConsumed[rname] += L;
    //};

    // ===============================================================================
    // 遍历 BAM 文件
    // ===============================================================================

    while(sam_read1(in, hdr, b) >= 0){
        const uint16_t f = b->core.flag;

        if(f & BAM_FSUPPLEMENTARY) continue;
        if(f & BAM_FSECONDARY) continue;
        if(f & BAM_FUNMAP) continue;

        std::string readName = bam_get_qname(b);
        uint16_t readFlag = f;

        if(b->core.tid >= 0){
            //consumeRefLen(hdr->target_name[b->core.tid], b->core.n_cigar, bam_get_cigar(b));
	    if(b->core.tid >= 0){
            int L = bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
            if(L > 0) {
                refConsumed[hdr->target_name[b->core.tid]] += L;
            }
        }
        }

        int main_mapq = b->core.qual;
        int main_nm = 0;
        if(uint8_t* nmTag = bam_aux_get(b, "NM"))
            main_nm = bam_aux2i(nmTag);

        int refLen1 = bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
        int matchLen1 = 0;
        uint32_t* c1 = bam_get_cigar(b);
        for(int i = 0; i < b->core.n_cigar; i++){
            int op = bam_cigar_op(c1[i]);
            int len = bam_cigar_oplen(c1[i]);
            if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
                matchLen1 += len;
        }

        if(!passMapqNm(main_mapq, main_nm)) continue;
        if(!passMatchFrac(matchLen1, std::max(1, refLen1))) continue;

        // ===============================================================================
        // 处理 split-read 证据（优先处理）
        // ===============================================================================

        uint8_t* saTag = bam_aux_get(b, "SA");
        bool hasSupplementEvidence = false;

        if(saTag && b->core.tid >= 0){
            std::string sa = bam_aux2Z(saTag);
            std::stringstream ss(sa);
            std::string item;

            std::string r1 = hdr->target_name[b->core.tid];
            int L1 = hdr->target_len[b->core.tid];
            int pos1 = b->core.pos + 1;
            bool rev1 = (b->core.flag & BAM_FREVERSE) != 0;
            ContigRegion reg1 = getContigRegion(pos1, L1);

            int readLen = getReadLength(b);

            std::string cigar1;
            {
                uint32_t* cig = bam_get_cigar(b);
                std::ostringstream oss;
                for(int i = 0; i < b->core.n_cigar; i++){
                    oss << bam_cigar_oplen(cig[i]) << bam_cigar_opchr(cig[i]);
                }
                cigar1 = oss.str();
            }
            ReadInterval interval1 = parseCigarReadInterval(cigar1, rev1, readLen);

            if(DEBUG_MODE){
                std::cerr << "\n=== Split-read: " << readName << " (len=" << readLen << ") ===\n";
                std::cerr << "Primary: " << r1 << " pos=" << pos1 << " rev=" << rev1
                          << " region=" << (reg1==START?"START":(reg1==END?"END":"MIDDLE"))
                          << " read[" << interval1.start << "-" << interval1.end << "]"
                          << " CIGAR=" << cigar1 << "\n";
            }

            while(std::getline(ss, item, ';')){
                if(item.empty()) continue;
                SAItem it = parseSAItem(item);
                if(!it.ok) continue;

                if(!passMapqNm(it.mapq, it.nm)) continue;

                int refLen2 = cigarRefLen(it.cigar);
                int matchLen2 = cigarMatchLen(it.cigar);
                if(!passMatchFrac(matchLen2, std::max(1, refLen2))) continue;

                std::string r2 = it.rname;

                if(r1 == r2) continue;

                auto tid_iter = name_to_tid.find(r2);
                if(tid_iter == name_to_tid.end()) continue;
                int tid2 = tid_iter->second;

                int L2 = hdr->target_len[tid2];
                int pos2 = it.pos;
                bool rev2 = it.isRev;
                ContigRegion reg2 = getContigRegion(pos2, L2);

                if(reg1 == MIDDLE || reg2 == MIDDLE) continue;

                ReadInterval interval2 = parseCigarReadInterval(it.cigar, rev2, readLen);

                if(DEBUG_MODE){
                    std::cerr << "SA: " << r2 << " pos=" << pos2 << " rev=" << rev2
                              << " region=" << (reg2==START?"START":(reg2==END?"END":"MIDDLE"))
                              << " read[" << interval2.start << "-" << interval2.end << "]"
                              << " CIGAR=" << it.cigar << "\n";
                }

                bool first1 = false;
                int maxGap = 150;
                int maxOverlap = 150;

                if(!canStitchReadIntervals(interval1, interval2, maxGap, maxOverlap, first1)) {
                    if(DEBUG_MODE) {
                        std::cerr << "  -> Cannot stitch: intervals too far apart or too much overlap\n";
                    }
                    continue;
                }

                if(DEBUG_MODE) {
                    std::cerr << "  -> Can stitch! " << (first1 ? "Primary first" : "SA first") << "\n";
                }

                bool found = false;
                std::string cL, cR;
                char oL_found, oR_found;

                for(char oL : {'+', '-'}) {
                    for(char oR : {'+', '-'}) {
                        if(checkSplitReadLayout(pos1, rev1, reg1, L1,
                                                pos2, rev2, reg2, L2,
                                                oL, oR, first1)) {
                            found = true;
                            cL = first1 ? r1 : r2;
                            cR = first1 ? r2 : r1;
                            oL_found = oL;
                            oR_found = oR;
                            goto split_found;
                        }
                    }
                }

                split_found:
                if(!found) {
                    if(DEBUG_MODE) {
                        std::cerr << "  -> No valid layout found\n";
                    }
                    continue;
                }

                if(DEBUG_MODE) {
                    std::cerr << "  -> Found valid layout: " << cL << "(" << oL_found
                              << ") -> " << cR << "(" << oR_found << ")\n";
                }

                Evidence ev;

                if(cL <= cR) {
                    ev.A = cL; ev.B = cR;
                    if(first1) {
                        ev.LA = L1; ev.LB = L2;
                        ev.posA = pos1; ev.posB = pos2;
                        ev.readOnAisRev = rev1;
                        ev.readOnBisRev = rev2;
                        ev.regA = reg1; ev.regB = reg2;
                        ev.mapqA = main_mapq; ev.nmA = main_nm;
                        ev.mapqB = it.mapq; ev.nmB = it.nm;
                    } else {
                        ev.LA = L2; ev.LB = L1;
                        ev.posA = pos2; ev.posB = pos1;
                        ev.readOnAisRev = rev2;
                        ev.readOnBisRev = rev1;
                        ev.regA = reg2; ev.regB = reg1;
                        ev.mapqA = it.mapq; ev.nmA = it.nm;
                        ev.mapqB = main_mapq; ev.nmB = main_nm;
                    }
                } else {
                    ev.A = cR; ev.B = cL;
                    if(first1) {
                        ev.LA = L2; ev.LB = L1;
                        ev.posA = pos2; ev.posB = pos1;
                        ev.readOnAisRev = rev2;
                        ev.readOnBisRev = rev1;
                        ev.regA = reg2; ev.regB = reg1;
                        ev.mapqA = it.mapq; ev.nmA = it.nm;
                        ev.mapqB = main_mapq; ev.nmB = main_nm;
                    } else {
                        ev.LA = L1; ev.LB = L2;
                        ev.posA = pos1; ev.posB = pos2;
                        ev.readOnAisRev = rev1;
                        ev.readOnBisRev = rev2;
                        ev.regA = reg1; ev.regB = reg2;
                        ev.mapqA = main_mapq; ev.nmA = main_nm;
                        ev.mapqB = it.mapq; ev.nmB = it.nm;
                    }
                }

                ev.readName = readName;
                ev.readFlag = readFlag;

                double score = 0.0;
                bool leftIsA = (ev.A == cL);
                char oL_eval = leftIsA ? oL_found : oR_found;
                char oR_eval = leftIsA ? oR_found : oL_found;

                if(computeLayoutScore(ev, leftIsA, oL_eval, oR_eval, score)) {
                    if(DEBUG_MODE) {
                        std::cerr << "  -> Passed eval with score=" << score << "\n";
                    }

                    LayoutKey key(cL, oL_found, cR, oR_found);
                    if(!OUTPUT_BOTH_ORDER && cR < cL) {
                        std::swap(cL, cR);
                        char newOL = (oR_found == '-') ? '+' : '-';
                        char newOR = (oL_found == '-') ? '+' : '-';
                        key = LayoutKey(cL, newOL, cR, newOR);
                    }

                    // 检查该连接是否在 fastg 预期中
                    ContigPair checkPair{cL, cR, oL_found, oR_found};
                    bool inFastg = (fastgPairs.find(checkPair) != fastgPairs.end());

                    auto &S = agg[key];
                    if(inFastg) {
                        S.supplementCount += 1;
                    } else {
                        S.supplementCountNoFastg += 1;
                    }
                    S.supportingReads.push_back(ReadInfo(ev.readName, ev.readFlag));

                    hasSupplementEvidence = true;
                }

                //refConsumed[r2] += std::max(0, refLen2);
            }
        }

        if(hasSupplementEvidence) {
            processedSupplementReads.insert(readName);
        }

        // ===============================================================================
        // 处理配对末端证据
        // ===============================================================================

        if(!hasSupplementEvidence && ENABLE_PAIRED && (f & BAM_FPAIRED) &&
           !(f & BAM_FMUNMAP) && b->core.mtid >= 0 && b->core.mtid != b->core.tid) {

            if(processedPairedReads.count(readName) > 0) {
                refConsumed[hdr->target_name[b->core.mtid]] += std::max(0, refLen1);
                continue;
            }

            std::string r1 = hdr->target_name[b->core.tid];
            std::string r2 = hdr->target_name[b->core.mtid];

            int L1 = hdr->target_len[b->core.tid];
            int L2 = hdr->target_len[b->core.mtid];

            int pos1 = b->core.pos + 1;
            int pos2 = b->core.mpos + 1;

            bool rev1 = (f & BAM_FREVERSE) != 0;
            bool rev2 = (f & BAM_FMREVERSE) != 0;

            ContigRegion reg1 = getContigRegion(pos1, L1);
            ContigRegion reg2 = getContigRegion(pos2, L2);

            if(reg1 != MIDDLE && reg2 != MIDDLE) {
                bool found = false;
                std::string cL, cR;
                char oL_found, oR_found;
                bool first1_found;

                for(int order = 0; order < 2; order++) {
                    for(char oL : {'+', '-'}) {
                        for(char oR : {'+', '-'}) {
                            bool first1 = (order == 0);

                            if(checkPairedEndLayout(pos1, rev1, reg1, L1,
                                                    pos2, rev2, reg2, L2,
                                                    oL, oR, first1)) {
                                found = true;
                                cL = first1 ? r1 : r2;
                                cR = first1 ? r2 : r1;
                                oL_found = oL;
                                oR_found = oR;
                                first1_found = first1;
                                goto paired_found;
                            }
                        }
                    }
                }

                paired_found:
                if(found) {
                    processedPairedReads.insert(readName);

                    Evidence ev;

                    if(cL <= cR) {
                        ev.A = cL; ev.B = cR;
                        if(first1_found) {
                            ev.LA = L1; ev.LB = L2;
                            ev.posA = pos1; ev.posB = pos2;
                            ev.readOnAisRev = rev1;
                            ev.readOnBisRev = rev2;
                            ev.regA = reg1; ev.regB = reg2;
                            ev.mapqA = main_mapq; ev.nmA = main_nm;
                            ev.mapqB = main_mapq; ev.nmB = main_nm;
                        } else {
                            ev.LA = L2; ev.LB = L1;
                            ev.posA = pos2; ev.posB = pos1;
                            ev.readOnAisRev = rev2;
                            ev.readOnBisRev = rev1;
                            ev.regA = reg2; ev.regB = reg1;
                            ev.mapqA = main_mapq; ev.nmA = main_nm;
                            ev.mapqB = main_mapq; ev.nmB = main_nm;
                        }
                    } else {
                        ev.A = cR; ev.B = cL;
                        if(first1_found) {
                            ev.LA = L2; ev.LB = L1;
                            ev.posA = pos2; ev.posB = pos1;
                            ev.readOnAisRev = rev2;
                            ev.readOnBisRev = rev1;
                            ev.regA = reg2; ev.regB = reg1;
                            ev.mapqA = main_mapq; ev.nmA = main_nm;
                            ev.mapqB = main_mapq; ev.nmB = main_nm;
                        } else {
                            ev.LA = L1; ev.LB = L2;
                            ev.posA = pos1; ev.posB = pos2;
                            ev.readOnAisRev = rev1;
                            ev.readOnBisRev = rev2;
                            ev.regA = reg1; ev.regB = reg2;
                            ev.mapqA = main_mapq; ev.nmA = main_nm;
                            ev.mapqB = main_mapq; ev.nmB = main_nm;
                        }
                    }

                    ev.readName = readName;
                    ev.readFlag = readFlag;

                    double score = 0.0;
                    bool leftIsA = (ev.A == cL);
                    char oL_eval = leftIsA ? oL_found : oR_found;
                    char oR_eval = leftIsA ? oR_found : oL_found;

                    if(computeLayoutScore(ev, leftIsA, oL_eval, oR_eval, score)) {
                        LayoutKey key(cL, oL_found, cR, oR_found);
                        if(!OUTPUT_BOTH_ORDER && cR < cL) {
                            std::swap(cL, cR);
                            char newOL = (oR_found == '-') ? '+' : '-';
                            char newOR = (oL_found == '-') ? '+' : '-';
                            key = LayoutKey(cL, newOL, cR, newOR);
                        }

                        // 检查该连接是否在 fastg 预期中
                        ContigPair checkPair{cL, cR, oL_found, oR_found};
                        bool inFastg = (fastgPairs.find(checkPair) != fastgPairs.end());

                        auto &S = agg[key];
                        if(inFastg) {
                            S.spanCount += 1;
                        } else {
                            S.spanCountNoFastg += 1;
                        }
                        S.supportingReads.push_back(ReadInfo(ev.readName, ev.readFlag));
                    }
                }
            }

            //refConsumed[r2] += std::max(0, refLen1);
        }
    }

    // ===============================================================================
    // 深度计算
    // ===============================================================================

    std::map<std::string, std::pair<double, int>> seg;
    for(int i = 0; i < hdr->n_targets; i++){
        std::string name = hdr->target_name[i];
        int L = hdr->target_len[i];
        if(L <= 0) continue;

        double consumed = 0.0;
        auto it = refConsumed.find(name);
        if(it != refConsumed.end()) consumed = it->second;

        double depth = consumed / std::max(1, L);
        double cnF = (avgDepth > 0.0) ? (depth / avgDepth) : 0.0;
        int cn = (int)std::floor(cnF + 0.5);

        seg[name] = {depth, cn};
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(in);

    // ===============================================================================
    // 输出结果
    // ===============================================================================

    std::ofstream out(outPath);
    if(!out){
        std::cerr << "Failed to open output " << outPath << "\n";
        return 1;
    }

    for(auto &kv : seg){
        out << "SEG " << kv.first << " " << kv.second.first << " " << kv.second.second << "\n";
    }

    for(const auto& kv : agg){
        const auto &lk = kv.first;
        const auto &S = kv.second;

        // 至少要有一种证据
        if(S.supplementCount == 0 && S.spanCount == 0 &&
           S.supplementCountNoFastg == 0 && S.spanCountNoFastg == 0) continue;

        // 应用最小计数过滤（总数）
        int totalCount = S.supplementCount + S.spanCount +
                         S.supplementCountNoFastg + S.spanCountNoFastg;
        if(totalCount < MIN_COUNT) continue;

        // 输出格式: JUNC left oL right oR supplementCount spanCount
        // 这里输出两个数字：支持 fastg 的计数 和 不支持 fastg 的计数
        out << "JUNC " << lk.left << " " << lk.oL << " "
            << lk.right << " " << lk.oR << " "
            << (S.supplementCount + S.spanCount + S.supplementCountNoFastg) << " "
            << (S.spanCountNoFastg);

        // 只在 DEBUG 模式下输出 reads 列表
        if(DEBUG_MODE) {
            out << " READS:";
            for(const auto& rinfo : S.supportingReads){
                out << " " << rinfo.name << "(" << rinfo.flag << ")";
            }
        }

        out << "\n";
    }

    return 0;
}
