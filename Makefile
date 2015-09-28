.PHONY: help \
        install-matlab install-octave install \
        uninstall-matlab uninstall-octave uninstall \
        test-matlab test-octave test

MATLAB?=matlab
OCTAVE?=octave

ROOTDIR=$(CURDIR)/surfing

ADDPATH="cd('$(ROOTDIR)');surfing_set_path();"
RMPATH="warning('off'); \
        rmpath(genpath('$(ROOTDIR)'));"
SAVEPATH="savepath();exit(0)"

BUILD="cd('$(CURDIR)'); \
       surfing_compile_mex;exit(0)"
INSTALL=$(ADDPATH)";"$(SAVEPATH)
UNINSTALL=$(RMPATH)";"$(SAVEPATH)


TEST=$(ADDPATH)";success=moxunit_runtests('${TESTDIR}');exit(~success);"

help:
	@echo "Usage: make <target>, where <target> is one of:"
	@echo "------------------------------------------------------------------"
	@echo "  install            to add surfing to the Matlab and GNU Octave"
	@echo "                     search paths, using whichever is present"
	@echo "  uninstall          to remove surfing from the Matlab and GNU"
	@echo "                     Octave search paths, using whichever is"
	@echo "                     present"
	@echo "  build              to compile mex code from the Matlab and GNU"
	@echo "                     Octave search paths, using whichever is"
	@echo "                     present"
	@echo ""
	@echo "  install-matlab     to add surfing to the Matlab search path"
	@echo "  install-octave     to add surfing to the GNU Octave search path"
	@echo "  uninstall-matlab   to remove surfing from the Matlab search path"
	@echo "  uninstall-octave   to remove surfing from the GNU Octave search"
	@echo "                     path"
	@echo "  build-matlab       to compile mex code using Matlab"
	@echo "  build-octave        to compile mex code using GNU Octave"
	@echo ""



MATLAB_BIN=$(shell which $(MATLAB))
OCTAVE_BIN=$(shell which $(OCTAVE))

ifeq ($(MATLAB_BIN),)
	# for Apple OSX, try to locate Matlab elsewhere if not found
    MATLAB_BIN=$(shell ls /Applications/MATLAB_R20*/bin/${MATLAB} 2>/dev/null | tail -1)
endif

MATLAB_RUN=$(MATLAB_BIN) -nojvm -nodisplay -nosplash -r
OCTAVE_RUN=$(OCTAVE_BIN) --no-gui --quiet --eval

build-matlab:
	@if [ -n "$(MATLAB_BIN)" ]; then \
		$(MATLAB_RUN) $(BUILD); \
	else \
		echo "matlab binary could not be found, skipping"; \
	fi;

build-octave:
	@if [ -n "$(OCTAVE_BIN)" ]; then \
		$(OCTAVE_RUN) $(BUILD); \
	else \
		echo "octave binary could not be found, skipping"; \
	fi;

build:
	@if [ -z "$(MATLAB_BIN)$(OCTAVE_BIN)" ]; then \
		@echo "Neither matlab binary nor octave binary could be found" \
		exit 1; \
	fi;
	$(MAKE) build-matlab
	$(MAKE) build-octave


install-matlab: build-matlab
	@if [ -n "$(MATLAB_BIN)" ]; then \
		$(MATLAB_RUN) $(INSTALL); \
	else \
		echo "matlab binary could not be found, skipping"; \
	fi;

install-octave: build-octave
	@if [ -n "$(OCTAVE_BIN)" ]; then \
		$(OCTAVE_RUN) $(INSTALL); \
	else \
		echo "octave binary could not be found, skipping"; \
	fi;

install:
	@if [ -z "$(MATLAB_BIN)$(OCTAVE_BIN)" ]; then \
		@echo "Neither matlab binary nor octave binary could be found" \
		exit 1; \
	fi;
	$(MAKE) install-matlab
	$(MAKE) install-octave


uninstall-matlab:
	@if [ -n "$(MATLAB_BIN)" ]; then \
		$(MATLAB_RUN) $(UNINSTALL); \
	else \
		echo "matlab binary could not be found, skipping"; \
	fi;

uninstall-octave:
	@if [ -n "$(OCTAVE_BIN)" ]; then \
		$(OCTAVE_RUN) $(UNINSTALL); \
	else \
		echo "octave binary could not be found, skipping"; \
	fi;

uninstall:
	@if [ -z "$(MATLAB_BIN)$(OCTAVE_BIN)" ]; then \
		@echo "Neither matlab binary nor octave binary could be found" \
		exit 1; \
	fi;
	$(MAKE) uninstall-matlab
	$(MAKE) uninstall-octave


test-matlab:
	@if [ -n "$(MATLAB_BIN)" ]; then \
		$(MATLAB_RUN) $(TEST); \
	else \
		echo "matlab binary could not be found, skipping"; \
	fi;

test-octave:
	@if [ -n "$(OCTAVE_BIN)" ]; then \
		$(OCTAVE_RUN) $(TEST); \
	else \
		echo "octave binary could not be found, skipping"; \
	fi;

test:
	@if [ -z "$(MATLAB_BIN)$(OCTAVE_BIN)" ]; then \
		@echo "Neither matlab binary nor octave binary could be found" \
		exit 1; \
	fi;
	$(MAKE) test-matlab
	$(MAKE) test-octave




