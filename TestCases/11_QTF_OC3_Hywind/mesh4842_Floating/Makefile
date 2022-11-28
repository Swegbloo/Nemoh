bindir=../../../bin/

.PHONY: run clean

# Build and run Sphere reference test case.
run:
			$(bindir)/preProc .
			$(bindir)/hydrosCal .
			$(bindir)/solver .
			$(bindir)/postProc .
			$(bindir)/QTFpreProc .
			$(bindir)/QTFsolver .
			$(bindir)/QTFpostProc .

# Delete results of the computations
clean:
			rm -rf results/
			rm -rf Normalvelocities.dat
			rm -rf /QTFPreprocOut
