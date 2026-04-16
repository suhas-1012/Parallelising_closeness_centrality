CXX = g++
CFLAGS = -O2 -std=c++17
OMPFLAGS = -fopenmp

SEQ_TARGETS = 01_seq_naive_bfs 02_seq_multisource_bfs 03_seq_bcc_reduced
PAR_TARGETS = 04_par_simple_mimd 05_par_mimd_msbfs 06_par_level_sync
DYN_TARGETS = 07_dynamic_shukla
NOVEL_TARGETS = 08_novel_bcc_spmm 09_novel_vdbcc

all: $(SEQ_TARGETS) $(PAR_TARGETS) $(DYN_TARGETS) $(NOVEL_TARGETS)

01_seq_naive_bfs: 01_seq_naive_bfs.cpp
	$(CXX) $(CFLAGS) -o $@ $<

02_seq_multisource_bfs: 02_seq_multisource_bfs.cpp
	$(CXX) $(CFLAGS) -o $@ $<

03_seq_bcc_reduced: 03_seq_bcc_reduced.cpp
	$(CXX) $(CFLAGS) -o $@ $<

04_par_simple_mimd: 04_par_simple_mimd.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

05_par_mimd_msbfs: 05_par_mimd_msbfs.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

06_par_level_sync: 06_par_level_sync.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

07_dynamic_shukla: 07_dynamic_shukla.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

08_novel_bcc_spmm: 08_novel_bcc_spmm.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

09_novel_vdbcc: 09_novel_vdbcc.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<		
	
clean:
	rm -f $(SEQ_TARGETS) $(PAR_TARGETS) $(DYN_TARGETS) $(NOVEL_TARGETS)

run_all: all
	@echo "========================================="
	@echo "01: Sequential Naive BFS"
	@echo "========================================="
	@./01_seq_naive_bfs
	@echo ""
	@echo "========================================="
	@echo "02: Sequential Multi-Source BFS"
	@echo "========================================="
	@./02_seq_multisource_bfs
	@echo ""
	@echo "========================================="
	@echo "03: Sequential BCC + R3/R4 Reduction"
	@echo "========================================="
	@./03_seq_bcc_reduced
	@echo ""
	@echo "========================================="
	@echo "04: Parallel Simple MIMD"
	@echo "========================================="
	@./04_par_simple_mimd
	@echo ""
	@echo "========================================="
	@echo "05: Parallel MIMD + MS-BFS"
	@echo "========================================="
	@./05_par_mimd_msbfs
	@echo ""
	@echo "========================================="
	@echo "06: Parallel Level-Sync Pull-Based"
	@echo "========================================="
	@./06_par_level_sync
	@echo ""
	@echo "========================================="
	@echo "07: Dynamic (Shukla-style, MIMD)"
	@echo "========================================="
	@./07_dynamic_shukla
	@echo ""
	@echo "========================================="
	@echo "08: Novel BCC-SpMM Hybrid (Dynamic)"
	@echo "========================================="
	@./08_novel_bcc_spmm
	@echo ""	@echo "========================================="
	@echo "09: Novel Vectorized Dynamic BCC Closeness Centrality"
	@echo "========================================="
	@./09_novel_vdbcc		

.PHONY: all clean run_all
