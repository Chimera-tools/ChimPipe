SHELL := /bin/bash

###################
#   Dependencies  #
###################
BEDTOOLS = bin/bedtools2-2.20.1
SAMTOOLS = bin/samtools-0.1.19
OVERLAP = bin/overlap-3.1
GEMTOOLS = bin/gemtools-1.7.1-i3/

###################
#    Makefiles    #
###################
BEDTOOLSMAKE = $(BEDTOOLS)/Makefile
SAMTOOLSMAKE = $(SAMTOOLS)/Makefile
OVERLAPMAKE = $(OVERLAP)/Makefile
GEMTOOLSMAKE = $(GEMTOOLS)/Makefile
MAKEFILES = $(BEDTOOLSMAKE) $(SAMTOOLSMAKE) $(OVERLAPMAKE) $(GEMTOOLSMAKE) 
#MAKEFILES = $(BEDTOOLSMAKE) $(SAMTOOLSMAKE) $(OVERLAPMAKE)  

#######################
#        Rules        #
#######################

all: make write run clean

compiled: write run clean
		
make: $(MAKEFILES)
	@echo
	@echo "Execute the Makefiles for the programs ChimPipe depends on"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@cd $(BEDTOOLS); make
	@cd $(SAMTOOLS); make
	@cd $(OVERLAP); make
	@cd $(GEMTOOLS); make
	@echo "done"
	
write: ChimPipe.sh
	@echo
	@echo "Write a temporary bash script"  
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "#!bin/bash" > temp.sh
	@echo "root=\$$(cd \$$(dirname \$$0); pwd -P )" >> temp.sh
	@echo "sed -i \"s|root=.*|root=\$$root|\" ChimPipe.sh" >> temp.sh
	@echo "done"
	
	
run:
	@echo
	@echo "Execute the script to add code to ChimPipe.sh to define an environmental" 
	@echo "variable with the directory where the pipeline is placed " 
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"	
	@$(shell bash temp.sh)
	@echo "done"
	
# Clean up
clean: temp.sh
	@echo
	@echo "4) Remove temporary bash script"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@rm -f temp.sh 
	@echo "done"


