SHELL := /bin/bash

###################
#   Dependencies  #
###################
BEDTOOLS = bin/bedtools2-2.20.1
SAMTOOLS = bin/samtools-0.1.19
OVERLAP = bin/Overlap-3.1
#GEMTOOLS = bin/gemtools-1.7.1

###################
#    Makefiles    #
###################
BEDTOOLSMAKE = $(BEDTOOLS)/Makefile
SAMTOOLSMAKE = $(SAMTOOLS)/Makefile
OVERLAPMAKE = $(OVERLAP)/Makefile
#GEMTOOLSMAKE = $(GEMTOOLS)/Makefile
MAKEFILESMAKE = $(BEDTOOLSMAKE) $(SAMTOOLSMAKE) $(OVERLAPMAKE) $(GEMTOOLSMAKE) 

#######################
#        Rules        #
#######################

all: make bash clean
				
make: $(MAKEFILES)
	@echo
	@echo "1) Execute the Makefiles for the programs ChimPipe depends on"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@cd $(BEDTOOLS); make
	@cd $(SAMTOOLS); make
	@cd $(OVERLAP); make
#	@cd $(GEMTOOLS); make
	
	
bash: ChimPipe.sh
	@echo
	@echo "2) Add code to ChimPipe.sh to define an environmental variable"
	@echo "with the directory where the pipeline is placed " 
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "create a temporal bash script"
	@echo "#!bin/bash" > temp.sh
	@echo "rootDir=\$$(cd \$$(dirname \$$0); pwd -P )" >> temp.sh
	@echo "sed -i \"165irootDir=\$$rootDir\" ChimPipe.sh" >> temp.sh
	@echo "execute it to add the code to ChimPipe"
	@$(shell bash ./temp.sh $<)
		
# Clean up
clean: temp.sh
	@echo
	@echo "3) Remove temporal bash script"
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
#	@rm -f temp.sh 
	@echo "finish"


