# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall

# Targets
TARGETS = HW_B HW_C HW_D HW_E HW_F
PYTHON_SCRIPTS = plot.py plot_F.py

# Output files
OUTPUT_FILES = Newton_output.txt Chebyshev_output.txt heart_10.txt heart_40.txt heart_160.txt

# Default target to build and run everything
all: $(TARGETS) run

# Compile each C++ program
HW_B: HW2_B.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

HW_C: HW2_C.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

HW_D: HW2_D.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

HW_E: HW2_E.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

HW_F: HW2_F.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# Run the compiled programs and Python scripts
run: $(TARGETS)
	@echo "Running HW_B"
	./HW_B
	@echo "Running HW_C"
	./HW_C
	@echo "Running plot.py"
	python3 plot.py
	@echo "Running HW_D"
	./HW_D
	@echo "Running HW_E"
	./HW_E
	@echo "Running HW_F"
	./HW_F
	@echo "Running plot_F.py"
	python3 plot_F.py

# Clean up all generated files
clean:
	rm -f $(TARGETS) $(OUTPUT_FILES)
	rm -f *.txt
