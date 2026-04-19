CXX      = g++
CFLAGS   = -O2 -std=c++17
OMPFLAGS = -fopenmp

SEQ_TARGETS   = 01_seq_naive_bfs 02_seq_multisource_bfs 03_seq_bcc_reduced
PAR_TARGETS   = 04_par_simple_mimd 05_par_mimd_msbfs 06_par_level_sync
DYN_TARGETS   = 07_dynamic_shukla
NOVEL_TARGETS = 08_novel_bcc_spmm 09_novel_vdbcc
GEN_TARGETS   = connected_graph_gen unconnected_graph_gen \
                bi_connected_graph_gen1 bi_connected_graph_gen2

ALL_TARGETS = $(SEQ_TARGETS) $(PAR_TARGETS) $(DYN_TARGETS) \
              $(NOVEL_TARGETS) $(GEN_TARGETS)

# ─── Default target ────────────────────────────────────────────────────────────
all: a $(ALL_TARGETS)

a:
	mkdir -p a

# ─── Sequential algorithms ─────────────────────────────────────────────────────
01_seq_naive_bfs: 01_seq_naive_bfs.cpp
	$(CXX) $(CFLAGS) -o $@ $<

02_seq_multisource_bfs: 02_seq_multisource_bfs.cpp
	$(CXX) $(CFLAGS) -o $@ $<

03_seq_bcc_reduced: 03_seq_bcc_reduced.cpp
	$(CXX) $(CFLAGS) -o $@ $<

# ─── Parallel algorithms ───────────────────────────────────────────────────────
04_par_simple_mimd: 04_par_simple_mimd.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

05_par_mimd_msbfs: 05_par_mimd_msbfs.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

06_par_level_sync: 06_par_level_sync.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

# ─── Dynamic algorithm ─────────────────────────────────────────────────────────
07_dynamic_shukla: 07_dynamic_shukla.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

# ─── Novel algorithms ──────────────────────────────────────────────────────────
08_novel_bcc_spmm: 08_novel_bcc_spmm.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

09_novel_vdbcc: 09_novel_vdbcc.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -o $@ $<

# ─── Graph generators ──────────────────────────────────────────────────────────
# NOTE: Rename bi_graph_gen1.cpp → bi_connected_graph_gen1.cpp
#             bi_graph_gen2.cpp → bi_connected_graph_gen2.cpp
#             Connected_graph_gen.cpp   → connected_graph_gen.cpp
#             Unconnected_graph_gen.cpp → unconnected_graph_gen.cpp

connected_graph_gen: connected_graph_gen.cpp
	$(CXX) $(CFLAGS) -o $@ $<

unconnected_graph_gen: unconnected_graph_gen.cpp
	$(CXX) $(CFLAGS) -o $@ $<

bi_connected_graph_gen1: bi_connected_graph_gen1.cpp
	$(CXX) $(CFLAGS) -o $@ $<

bi_connected_graph_gen2: bi_connected_graph_gen2.cpp
	$(CXX) $(CFLAGS) -o $@ $<

# ─── Utility targets ───────────────────────────────────────────────────────────
clean:
	rm -f $(ALL_TARGETS)

clean_data:
	rm -f experiment.csv connected_graph.csv unconnected_graph.csv \
	      biconnected_graph.csv a/*.csv

clean_all: clean clean_data

seq:  a $(SEQ_TARGETS)
par:  a $(PAR_TARGETS)
dyn:  a $(DYN_TARGETS)
novel: a $(NOVEL_TARGETS)
gen:  $(GEN_TARGETS)

.PHONY: all clean clean_data clean_all seq par dyn novel gen