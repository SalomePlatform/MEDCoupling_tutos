TO = beamer
# pdf html js

exec = pandoc
opt_common = -s -N # --embed-resources

ifeq ($(TO),beamer)
format = pdf
opt_beamer = --slide-level=2 -t beamer -V theme=CEA23
exec_opt = $(exec) $(opt_common) $(filename) $(opt_beamer)
prefix = pres_

else ifeq ($(TO),pdf)
format = pdf
exec_opt = $(exec) $(opt_common) $(filename)
prefix =

endif

../install/$(prefix)%.$(format): ../src/%.md
	# creates dir if not already created
	mkdir -p ../install
	$(exec_opt) -o $@ $<

.PHONY: pres_1
pres_1: ../install/$(prefix)1_medcoupling.$(format)

.PHONY: pres_2
pres_2: ../install/$(prefix)2_medfile_io.$(format)

.PHONY: pres_3
pres_3: ../install/$(prefix)3_remapper.$(format)

.PHONY: pres_4
pres_4: ../install/$(prefix)4_parallel.$(format)

.PHONY: pres_resume
pres_resume: ../install/$(prefix)resume_francais.$(format)

.PHONY: all
all: pres_1 pres_2 pres_3 pres_4

clean:
	-rm ../install/*pdf
	-rm ../install/*html
	-rm *.err

clean_all: clean
	cp Makefile /tmp/
	-rm -r *
	cp /tmp/Makefile .
