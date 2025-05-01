//
// Created by gzpan2 on 4/25/2024.
//
#include <htslib/sam.h>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath> // For floor function/min
const int MAX_REVERSE = 150;
const int MAX_END = 300;
struct ContigPair {
    std::string ref1;
    std::string ref2;
    int orientation; // 0: ++, 1: +-, 2: -+, 3: --
    bool operator<(const ContigPair& other) const {
        if (ref1 != other.ref1) return ref1 < other.ref1;
        if (ref2 != other.ref2) return ref2 < other.ref2;
        return orientation < other.orientation;
    }
};
struct SegmentedContig {
    std::string name;
    double depth;
    int copyNum;
};
std::string orientationToString(int orientation) {
    switch (orientation) {
        case 0: return "++";
        case 1: return "+-";
        case 2: return "-+";
        case 3: return "--";
        default: return "??";
    }
}

int getEdgeLen(std::string edge) {
    std::istringstream iss(edge);
    std::string token;
    std::vector<std::string> tokens;

    // Split the string by underscores
    while (std::getline(iss, token, '_')) {
        tokens.push_back(token);
    }

    // Check if there are at least 4 elements
//    if (tokens.size() >= 4) {
//        std::cout << "The 4th value is: " << tokens[3] << std::endl; // Index 3 because indexing starts at 0
//    } else {
//        std::cout << "There are less than 4 values in the string." << std::endl;
//    }
    return std::stoi(tokens[3]);
//    return -1;
}
std::map<ContigPair, int> parseFastgFile(const std::string& filename) {
    std::ifstream inFile(filename);
    std::string line;
    std::map<ContigPair, int> contigPairs;

    while (getline(inFile, line)) {
//        if (line.front() == '>') {  // FastG headers start with '>'
        std::stringstream ss(line);  // Skip the '>'
        std::string fullName, contigName, linkedContig;
        std::getline(ss, fullName, ';');
        std::stringstream ss2(fullName);
        std::getline(ss2, contigName, ':');  // Get the contig name before ':'

        // Now read the relationships
        while (std::getline(ss2, linkedContig, ',')) {
            // Assume format is like "contig2+:contig3-"
            int orientation = 0;
            if(contigName.back() != '\''){
                if (linkedContig.back() != '\''){
                    orientation = 0;
                } else {
                    linkedContig.pop_back();
                    orientation = 1;
                }
            } else {
                contigName.pop_back();
                if (linkedContig.back() != '\''){
                    orientation = 2;
                } else {
                    linkedContig.pop_back();
                    orientation = 3;
                }
            }
//
//                    if (linkedContig.back() != '\'') {
//                    orientation = (linkedContig[linkedContig.size() - 2] == '\'') ? 0 : 2;
//                } else {
//                    orientation = (linkedContig[linkedContig.size() - 2] == '\'') ? 1 : 3;
//                }

//                std::string ref2 = linkedContig.substr(0, linkedContig.size() - 2);
            ContigPair pair{contigName, linkedContig, orientation};
            contigPairs[pair] = 0;  // Initialize or update the count/metric
        }
//        }
    }

    return contigPairs;
}
void get_values(int supp, int fai, int r1r2, std::array<int, 7>& result) {
    if(supp && fai && r1r2) {
        result[0]++;
    } else if (supp && fai) {
        result[1]++;
    } else if (fai && r1r2) {
        result[2]++;
    } else if (supp && r1r2) {
        result[3]++;
    } else if (supp) {
        result[4]++;
    } else if (fai) {
        result[5]++;
    } else if (r1r2) {
        result[6]++;
    }
}
int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <BAM file> <assembly.fastg.fai> <Out file> <Average Depth>" << std::endl;
        return 1;
    }

    samFile *in = sam_open(argv[1], "r");
    std::map<ContigPair, int> fastgPairs = parseFastgFile(argv[2]);
    std::ofstream outFile(argv[3]);
    auto averageDepth = std::stod(argv[4]);

    if (in == nullptr) {
        std::cerr << "Failed to open BAM file " << argv[1] << std::endl;
        return 1;
    }

    bam_hdr_t *hdr = sam_hdr_read(in);
    if (hdr == nullptr) {
        std::cerr << "Failed to read header from BAM file" << std::endl;
        sam_close(in);
        return 1;
    }
    std::string tags;

    bam1_t *aln = bam_init1();
    std::map<ContigPair, std::array<int, 7>> contigPairs;
    std::map<std::string, double> totalBaseCounts;
    std::map<std::string, SegmentedContig> segmentedContigs;
    while (sam_read1(in, hdr, aln) >= 0) {
        if (aln->core.tid >= 0) {  // Ensure we have a valid target ID
            const char* contigName = hdr->target_name[aln->core.tid];
            int alignedLength = bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln));
            totalBaseCounts[contigName] += alignedLength;
        }
        if (strcmp(bam_get_qname(aln), "virulent:NC_020855.1_0-33574") == 0) {
            int tmp = 0;
        }
        int supp = 0;
        int fai = 0;
        int r1r2 = 0;
        bool is_continue = false;
        // Inside the loop where you read alignments
        uint8_t *sa_tag = bam_aux_get(aln, "SA");
        uint8_t *xa_tag = bam_aux_get(aln, "XA");
        if(xa_tag) continue;
        if (sa_tag) {
            // Handle supplementary alignments (split reads)
            uint32_t read_refID = aln->core.tid;
            std::string read_refName = hdr->target_name[read_refID];
            bool read_reversed = bam_is_rev(aln);
//            uint8_t *sa_tag = bam_aux_get(aln, "SA");
//            if (sa_tag) {
            std::string sa_str = bam_aux2Z(sa_tag);
            std::stringstream ss(sa_str);
            std::string item;
            while (std::getline(ss, item, ';')) {
                std::stringstream item_stream(item);
                std::string sa_refName, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm;
                std::getline(item_stream, sa_refName, ',');
                std::getline(item_stream, sa_pos, ',');
                std::getline(item_stream, sa_strand, ',');
                std::getline(item_stream, sa_cigar, ',');
                std::getline(item_stream, sa_mapq, ',');
                std::getline(item_stream, sa_nm, ',');

                if (!sa_refName.empty() && sa_refName != read_refName) {
                    bool next_reversed = sa_strand == "-";
                    // Calculate orientation based on strand
                    int orientation = 0;
                    if (read_reversed == next_reversed) {
                        if(std::stoi(sa_pos) <= 10) {
                            orientation = 0;
                        } else {
                            orientation = 3;
                        }
                    } else {
                        orientation = 2;
                    }
//                        int orientation = (read_reversed ? 2 : 0) + (sa_strand == "-" ? 1 : 0);
                    ContigPair pair = {read_refName, sa_refName, orientation};
//                    contigPairs[pair]++;
                    supp = 1;
//                    tags.append(" supp");
                    if (aln->core.mtid != -1) {
                        std::string refName2 = hdr->target_name[aln->core.mtid];
                        if (refName2 == sa_refName) {
//                            tags.append(" pair");
                            r1r2 = 1;
                        } else {
                            continue;
                        }
                    }
                    if (fastgPairs.count(pair) != 0) {
//                        tags.append(" fastg");
                        fai = 1;
                    }
                    if (contigPairs.count(pair) != 0) {
                        get_values(supp, fai, r1r2, contigPairs[pair]);
                    } else {
                        std::array<int, 7> array{};
                        get_values(supp, fai, r1r2, array);
                        contigPairs[pair] = array;
                    }
                } else {
                    is_continue = true;
                }
            }
//            }
        }
        if ((!sa_tag || is_continue) && !(aln->core.flag & BAM_FUNMAP) && aln->core.mtid != -1) {
            auto m = bam_get_qname(aln);
            if (strcmp(bam_get_qname(aln), "virulent:NC_020855.1_0-33574") == 0) {
                int tmp = 0;
            }
//            if (aln->core.tid == aln->core.mtid) {
////                uint32_t pos = aln->core.pos + 1; // Position of the alignment
////                uint32_t mpos = aln->core.mpos + 1; // Mate's position
////                std::string refName1 = hdr->target_name[aln->core.tid];
////                const char* qname = bam_get_qname(aln); // Query name of the read
////                // Check flags for orientation
////                bool is_read1_reverse = bam_is_rev(aln);
////                bool is_read2_reverse = bam_is_mrev(aln);
////                if (((pos < MAX_END && mpos > hdr->target_len[aln->core.tid] - MAX_END) &&
////                     !is_read1_reverse & is_read2_reverse) ||
////                    ((mpos < MAX_END && pos > hdr->target_len[aln->core.tid] - MAX_END) &&
////                     is_read1_reverse & !is_read2_reverse)) {
////                    int orientation = 0;
////                    ContigPair pair = {refName1, refName1, orientation};
////                    r1r2 = 1;
////                    if (fastgPairs.count(pair) != 0) {
////                        fai = 1;
////                    }
////                    if (contigPairs.count(pair) != 0) {
////                        get_values(supp, fai, r1r2, contigPairs[pair]);
////                    } else {
////                        std::array<int, 7> array{};
////                        get_values(supp, fai, r1r2, array);
////                        contigPairs[pair] = array;
////                    }
//////                std::cout << "Read pair (" << qname << ") spans the head and tail of contig "
//////                          << hdr->target_name[aln->core.tid] << " with proper orientation." << std::endl;
////                }
//            } else
            if (aln->core.flag & BAM_FREAD1) {
                // Handle properly paired reads on different contigs (span reads)
                std::string refName1 = hdr->target_name[aln->core.tid];
                std::string refName2 = hdr->target_name[aln->core.mtid];
                uint32_t pos = aln->core.pos + 1; // Position of the alignment
                uint32_t mpos = aln->core.mpos + 1; // Mate's position
//                ContigPair pair;
                if (refName2 == "EDGE_98062_length_57763_cov_76.403376" &&
                    refName1 == "EDGE_104482_length_333_cov_1232.618705") {
                    int tmp = 0;
                }
                bool read_reversed = bam_is_rev(aln);
                bool mate_reversed = bam_is_mrev(aln);
                int orientation = -1;
                int tid_pref = std::min(int(hdr->target_len[aln->core.tid] / 2), MAX_END);
                int tid_suff = std::max(int(hdr->target_len[aln->core.tid] / 2),
                                        int(hdr->target_len[aln->core.tid]) - MAX_END);

                int mtid_pref = std::min(int(hdr->target_len[aln->core.mtid] / 2), MAX_END);
                int mtid_suff = std::max(int(hdr->target_len[aln->core.mtid] / 2),
                                         int(hdr->target_len[aln->core.mtid]) - MAX_END);

                int pd = pos < tid_pref ? 0 : (pos > tid_suff ? 1 : -1);
                int mpd = (mpos < mtid_pref) ? 0 : (mpos > mtid_suff ? 1 : -1);
//                if ref or mref is two shot
                std::string pref1, pref2;
                if (pd == 1) {
                    if (!read_reversed) {
                        if (mpd == 0) {
                            if (mate_reversed) {
                                orientation = 0;
                            } else {
//                                if mreference is too short
                                if ((hdr->target_len[aln->core.mtid] < (mpos + aln->core.l_qseq)) || (hdr->target_len[aln->core.mtid] - (mpos + aln->core.l_qseq) < MAX_REVERSE)) {
                                    orientation = 1;
                                }
                            }
                        } else if (mpd == 1) {
                            if (!mate_reversed) {
                                orientation = 1;
                            } else {
//                                if mreference is too short
                                if (mpos < MAX_REVERSE) {
                                    orientation = 0;
                                }
                            }
                        }
//                        pref1 = refName1;
//                        pref2 = refName2;
                    } else if (pos < MAX_REVERSE) {
                        if (mpd == 0) {
                            if (mate_reversed) {
                                orientation = 2;
                            } else {
//                                if mreference is too short
                                if ((hdr->target_len[aln->core.mtid] < (mpos + aln->core.l_qseq)) || (hdr->target_len[aln->core.mtid] - (mpos + aln->core.l_qseq) < MAX_REVERSE)) {
                                    orientation = 3;
                                }
                            }
                        } else if (mpd == 1) {
                            if (!mate_reversed) {
                                orientation = 3;
                            } else {
//                                if mreference is too short
                                if (mpos < MAX_REVERSE) {
                                    orientation = 2;
                                }
                            }
                        }
//                        pref1 = refName1;
//                        pref2 = refName2;
                    }
                } else {
                    if (read_reversed) {
                        if (mpd == 0) {
                            if (mate_reversed) {
                                orientation = 2;
                            } else if ((hdr->target_len[aln->core.mtid] < (mpos + aln->core.l_qseq)) ||
                                       (hdr->target_len[aln->core.mtid] - (mpos + aln->core.l_qseq) < MAX_REVERSE)) {
                                orientation = 3;
                            }
                        } else if (mpd == 1) {
                            if (!mate_reversed) {
                                orientation = 3;
                            } else if (mpos < MAX_REVERSE) {
                                {
                                    orientation = 2;
                                }
                            }
                        } else if ((hdr->target_len[aln->core.tid] < (pos + aln->core.l_qseq)) ||
                                   hdr->target_len[aln->core.tid] - (pos + aln->core.l_qseq) < MAX_REVERSE) {
                            if (mpd == 0) {
                                if (mate_reversed) {
                                    orientation = 0;
                                } else if ((hdr->target_len[aln->core.mtid] < (mpos + aln->core.l_qseq)) ||
                                           (hdr->target_len[aln->core.mtid] - (mpos + aln->core.l_qseq) <
                                            MAX_REVERSE)) {
                                    orientation = 1;
                                }
                            } else if (mpd == 1) {
                                if (!mate_reversed) {
                                    orientation = 1;
                                } else if (mpos < MAX_REVERSE) {
                                    orientation = 0;
                                }
                            }
//                            pref1 = refName1;
//                            pref2 = refName2;
                        }
                    }
                }

                    if (orientation != -1) {
                        pref1 = refName1;
                        pref2 = refName2;
                        ContigPair pair = {pref1, pref2, orientation};
                        r1r2 = 1;
//            if(!xa_tag) {
//                tags.append(" no_xa");
//            }
                        if (fastgPairs.count(pair) != 0) {
//                tags.append(" fastg");
                            fai = 1;
                        }
                        if (contigPairs.count(pair) != 0) {
                            get_values(supp, fai, r1r2, contigPairs[pair]);
                        } else {
                            std::array<int, 7> array{};
                            get_values(supp, fai, r1r2, array);
                            contigPairs[pair] = array;
                        }
                    }

//                if ((mpos < MAX_END && pos > hdr->target_len[aln->core.tid] - MAX_END)) {
//                    if (!read_reversed && mate_reversed) {
//                        orientation = 0;
//                    } else if (!read_reversed && !mate_reversed) {
//                        orientation = 1;
//                    } else if (read_reversed && !mate_reversed) {
//                        orientation = 3;
//                    } else {
//                        orientation = 1;
//                    }
////            int orientation = (read_reversed ? 2 : 0) + (mate_reversed ? 1 : 0);
//                    ContigPair pair = {refName1, refName2, orientation};
//                    r1r2 = 1;
////            if(!xa_tag) {
////                tags.append(" no_xa");
////            }
//                    if (fastgPairs.count(pair) != 0) {
////                tags.append(" fastg");
//                        fai = 1;
//                    }
//                    if (contigPairs.count(pair) != 0) {
//                        get_values(supp, fai, r1r2, contigPairs[pair]);
//                    } else {
//                        std::array<int, 7> array{};
//                        get_values(supp, fai, r1r2, array);
//                        contigPairs[pair] = array;
//                    }
//                }
//                else if ((pos < MAX_END && mpos > hdr->target_len[aln->core.mtid] - MAX_END)) {
//                    if (!mate_reversed && read_reversed) {
//                        orientation = 0;
//                    } else if (!mate_reversed && !read_reversed) {
//                        orientation = 1;
//                    } else if (mate_reversed && !read_reversed) {
//                        orientation = 3;
//                    } else {
//                        orientation = 1;
//                    }
//                    ContigPair pair = {refName2, refName1, orientation};
//                    r1r2 = 1;
////            if(!xa_tag) {
////                tags.append(" no_xa");
////            }
//                    if (fastgPairs.count(pair) != 0) {
////                tags.append(" fastg");
//                        fai = 1;
//                    }
//                    if (contigPairs.count(pair) != 0) {
//                        get_values(supp, fai, r1r2, contigPairs[pair]);
//                    } else {
//                        std::array<int, 7> array{};
//                        get_values(supp, fai, r1r2, array);
//                        contigPairs[pair] = array;
//                    }
//                }

//            contigPairs[pair]++;
//            tags.append(" pair");


            }
        }
    }
    // Calculate depth and estimate copy number
    for (int i = 0; i < hdr->n_targets; ++i) {
        std::string contigName = hdr->target_name[i];
        int contigLength = hdr->target_len[i];
        double depth = totalBaseCounts[contigName] / contigLength;
        double copyNum = static_cast<double >(depth / averageDepth); // Example: assume 1 copy per 10x depth
        // Convert copyNum to int with conditional rounding
        int intCopyNum = static_cast<int>(floor(copyNum));
        if (copyNum - intCopyNum > 0.7) {
            intCopyNum += 1;
        }
        segmentedContigs[contigName] = {contigName, depth, intCopyNum};
    }

    bam_destroy1(aln);
    bam_hdr_destroy(hdr);
    sam_close(in);

    if (!outFile) {
        std::cerr << "Failed to open file for writing." << std::endl;
        return 1;
    }
    for (const auto& entry : segmentedContigs) {
        outFile << "SEG " << entry.second.name << " " << entry.second.depth << " " << entry.second.copyNum << std::endl;
    }

    // Write the contents previously going to std::cout to the file
//    for (auto &pair : contigPairs) {
//        outFile <<"JUNC "<< pair.first.ref1 << " " <<orientationToString(pair.first.orientation)[0] <<" "<< pair.first.ref2 <<" " <<orientationToString(pair.first.orientation)[1]
//                << " " << pair.second[0] <<" "<<pair.second[1]<<" "<<pair.second[2]<<" "<<pair.second[3]<<" "<<pair.second[4]<<" "<<pair.second[5]<<" "<<pair.second[6] << std::endl;
//    }
    for (auto &pair : contigPairs) {
//        if (pair.second[0] + pair.second[1] + pair.second[2] <= 5 && pair.second[6] + pair.second[4] <= 10) {
//            continue;
//        }
        if (pair.second[0] + pair.second[1] + pair.second[2] + pair.second[6] + pair.second[4]< 5) {
            continue;
        }
//        if (getEdgeLen(pair.first.ref1) < 150 || getEdgeLen(pair.first.ref2) < 150) {
//            continue;
//        }
//        if (pair.second[0] + pair.second[1] + pair.second[2] == 0 && (getEdgeLen(pair.first.ref1) < 300 || getEdgeLen(pair.first.ref2) < 300)) {
//            continue;
//        }
        outFile <<"JUNC "<< pair.first.ref1 << " " <<orientationToString(pair.first.orientation)[0] <<" "<< pair.first.ref2 <<" " <<orientationToString(pair.first.orientation)[1]
                << " " << pair.second[0] + pair.second[1] + pair.second[2] <<" "<< pair.second[6] + pair.second[4] << std::endl;
    }
    // for(auto &p:fastgPairs) {
    //     if (contigPairs.count(p.first) == 0) {
    //         outFile <<"JUNC "<< p.first.ref1 << " " <<orientationToString(p.first.orientation)[0] <<" "<< p.first.ref2 <<" " <<orientationToString(p.first.orientation)[1]
    //                 << " " << 0 <<" "<< -1 << std::endl;
    //     }
    // }

    // Close the file stream
    outFile.close();

    return 0;
}