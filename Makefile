SHELL := /bin/bash

#######################
#        Rules        #
#######################
all: write run clean

write:		
	@echo "#!bin/bash" > temp.sh
	@echo "root=\$$(cd \$$(dirname \$$0); pwd -P )" >> temp.sh
	@echo "sed -i \"s|root=.*|root=\$$root|\" ChimPipe.sh" >> temp.sh
		
run:
	@$(shell bash temp.sh)

# Clean up
clean: temp.sh
	@echo
	@rm -f temp.sh 
	@echo "Done!"
	@echo

