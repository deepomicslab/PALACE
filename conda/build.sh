#!/bin/bash
set -euxo pipefail

if [[ ! -d "bin" ]]; then
    echo "Error: Source directory 'bin' not found."
    exit 1
fi

echo "Installing files to ${PREFIX}..."

mkdir -p ${PREFIX}/bin
mkdir -p ${PREFIX}/share/palace/scripts

echo "Installing binaries..."
# 使用 Conda 传入的编译器环境变量 ($CXX, $CXXFLAGS, $LDFLAGS) 进行安全编译
${CXX} ${CXXFLAGS} -O2 -std=c++17 \
    -I"${PREFIX}/include" \
    -L"${PREFIX}/lib" \
    bin/generate_graph.cpp \
    -lhts -lz -lbz2 -llzma -lcurl -lpthread \
    -o bin/generateGraph ${LDFLAGS}

${CXX} ${CXXFLAGS} -O2 -std=c++17 \
    -I"${PREFIX}/include" \
    -L"${PREFIX}/lib" \
    bin/extract_ref.cpp \
    -lhts -lz -lbz2 -llzma -lcurl -lpthread \
    -o bin/eref ${LDFLAGS}

cp bin/generateGraph ${PREFIX}/bin/
cp bin/eref ${PREFIX}/bin/
cp bin/matching ${PREFIX}/bin/
cp palace ${PREFIX}/bin/
chmod +x ${PREFIX}/bin/*

# 将 Python Cython/C 扩展的编译提前到 build 阶段
echo "Building Python extensions..."
if [[ -d "share/palace/scripts" ]] && [[ -f "share/palace/scripts/setup.py" ]]; then
    pushd share/palace/scripts
    ${PYTHON} setup.py build_ext --inplace
    popd
fi

echo "Installing scripts..."
if [[ -d "share/palace/scripts" ]]; then
    find share/palace/scripts/ -maxdepth 1 -type f \( -name "*.py" -o -name "*.so" -o -name "*.sh" -o -name "*.pyx" \) \
        -exec cp {} ${PREFIX}/share/palace/scripts/ \;

    chmod +x ${PREFIX}/share/palace/scripts/*.py 2>/dev/null || true
    chmod +x ${PREFIX}/share/palace/scripts/*.sh 2>/dev/null || true
fi

echo "Installing config file..."
if [[ -f "config/config.txt" ]]; then
    cp config/config.txt ${PREFIX}/share/palace/
fi

# 将 palace-env 包装器的创建提前到 build 阶段
echo "Creating environment wrapper for palace..."
cat > "${PREFIX}/bin/palace-env" << 'EOF'
#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
export PATH="${SCRIPT_DIR}/../share/palace/scripts:${PATH}"
exec "${SCRIPT_DIR}/palace" "$@"
EOF
chmod +x "${PREFIX}/bin/palace-env"

echo "Build script finished successfully."