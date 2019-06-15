# The following two variable will be commandline determined by TA
# For testing, you could uncomment them.
SRIPATH ?= /Users/apple/Desktop/srilm
MACHINE_TYPE ?= macosx
LM ?= bigram.lm

CXX = g++
CXXFLAGS = -O3 -I$(SRIPATH)/include -w
vpath lib%.a $(SRIPATH)/lib/$(MACHINE_TYPE)

TARGET = mydisambig
SRC = mydisambig.cpp
OBJ = $(SRC:.cpp=.o)
TO = ZhuYin-Big5.map
FROM = Big5-ZhuYin.map
.PHONY: all clean map run

all: $(TARGET)

OUTPUT_PATH := result2


$(TARGET): $(OBJ) -loolm -ldstruct -lmisc
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
run:
	[ -d result2 ] || mkdir -p result2;
	for i in $(shell seq 1 10) ; do \
	    echo "Running $$i.txt"; \
	    ./mydisambig testdata/$$i.txt  > result2/$$i.txt 2>nul; \
	done;
map:
	echo "Convert Big5-ZhuYin.map to ZhuYin-Big5.map..."
	python3 mapping.py
clean:
	$(RM) $(OBJ) $(TARGET)
diff:
# 	if diff result1/1.txt result1/1.txt; then echo "same"; else echo "qqq"; fi

	for i in $(shell seq 1 10) ; do \
	    echo "Running $$i.txt"; \
	    diff result1/$$i.txt result2/$$i.txt || echo "difference detected"; \
	done;
