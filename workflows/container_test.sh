#!/bin/bash
# container_test.sh - Test bioinformatics containers for Day 7 training
# Usage: ./container_test.sh

echo "🧬 Testing Bioinformatics Containers for MTB Pipeline"
echo "=================================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to test container
test_container() {
    local container=$1
    local command=$2
    local tool_name=$3
    
    echo -e "${YELLOW}Testing ${tool_name}...${NC}"
    
    if singularity exec docker://${container} ${command} > /dev/null 2>&1; then
        echo -e "${GREEN}✅ ${tool_name} container works${NC}"
        return 0
    else
        echo -e "${RED}❌ ${tool_name} container failed${NC}"
        return 1
    fi
}

# Test FastQC
test_container "biocontainers/fastqc:v0.11.9" "fastqc --version" "FastQC"

# Test Trimmomatic
test_container "staphb/trimmomatic:0.39" "trimmomatic -version" "Trimmomatic"

# Test SPAdes
test_container "staphb/spades:3.15.4" "spades.py --version" "SPAdes"

# Test Prokka
test_container "staphb/prokka:1.14.6" "prokka --version" "Prokka"

# Test AMRFinderPlus
test_container "staphb/amrfinderplus:3.10.23" "amrfinder --version" "AMRFinderPlus"

# Test TB-Profiler
test_container "staphb/tbprofiler:4.1.1" "tb-profiler version" "TB-Profiler"

# Test QUAST
test_container "staphb/quast:5.0.2" "quast.py --version" "QUAST"

# Test MLST
test_container "staphb/mlst:2.22.0" "mlst --version" "MLST"

# Test MultiQC
test_container "ewels/multiqc:v1.12" "multiqc --version" "MultiQC"

echo ""
echo "🎉 Container testing complete!"
echo ""
echo "💡 Tips:"
echo "- If any containers failed, check your internet connection"
echo "- Singularity will download containers on first use"
echo "- Use 'singularity cache clean' to clear cached containers"
echo "- Run 'singularity exec docker://container:tag command --help' for tool help"
