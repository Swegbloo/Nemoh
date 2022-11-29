SHELL=/bin/bash

verdir=./

.PHONY: run clean
run: 
	@echo ""
	@echo "Sphere"
	@$(MAKE) --silent -C $(verdir)/1_Sphere/                     test
	@echo ""
	@echo "Sphere using y-symmetry"
	@$(MAKE) --silent -C $(verdir)/2_SymmetricSphere/            test
	@echo ""
	@echo "Sphere in finite depth"
	@$(MAKE) --silent -C $(verdir)/3_FiniteDepthSphere/          test
	@echo ""
	@echo "Sphere in finite depth using y-symmetry"
	@$(MAKE) --silent -C $(verdir)/4_SymmetricFiniteDepthSphere/ test
	@echo ""
	@echo "Alien sphere"
	@$(MAKE) --silent -C $(verdir)/5_AlienSphere/                test

clean:
	@$(MAKE) --silent -C $(verdir)/1_Sphere/                     clean
	@$(MAKE) --silent -C $(verdir)/2_SymmetricSphere/            clean
	@$(MAKE) --silent -C $(verdir)/3_FiniteDepthSphere/          clean
	@$(MAKE) --silent -C $(verdir)/4_SymmetricFiniteDepthSphere/ clean
	@$(MAKE) --silent -C $(verdir)/5_AlienSphere/                clean

