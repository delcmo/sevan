sevan_INC_DIRS := $(shell find $(SEVAN_DIR)/include -type d -not -path "*/.svn*")
sevan_INCLUDE  := $(foreach i, $(sevan_INC_DIRS), -I$(i))

libmesh_INCLUDE := $(sevan_INCLUDE) $(libmesh_INCLUDE)

sevan_LIB := $(SEVAN_DIR)/libsevan-$(METHOD).la

sevan_APP := $(SEVAN_DIR)/sevan-$(METHOD)

# source files
sevan_srcfiles    := $(shell find $(SEVAN_DIR)/src -name "*.C" -not -name main.C)
sevan_csrcfiles   := $(shell find $(SEVAN_DIR)/src -name "*.c")
sevan_fsrcfiles   := $(shell find $(SEVAN_DIR)/src -name "*.f")
sevan_f90srcfiles := $(shell find $(SEVAN_DIR)/src -name "*.f90")

# object files
sevan_objects := $(patsubst %.C, %.$(obj-suffix), $(sevan_srcfiles))
sevan_objects += $(patsubst %.c, %.$(obj-suffix), $(sevan_csrcfiles))
sevan_objects += $(patsubst %.f, %.$(obj-suffix), $(sevan_fsrcfiles))
sevan_objects += $(patsubst %.f90, %.$(obj-suffix), $(sevan_f90srcfiles))

# plugin files
sevan_plugfiles    := $(shell find $(SEVAN_DIR)/plugins/ -name "*.C" 2>/dev/null)
sevan_cplugfiles   := $(shell find $(SEVAN_DIR)/plugins/ -name "*.c" 2>/dev/null)
sevan_fplugfiles   := $(shell find $(SEVAN_DIR)/plugins/ -name "*.f" 2>/dev/null)
sevan_f90plugfiles := $(shell find $(SEVAN_DIR)/plugins/ -name "*.f90" 2>/dev/null)

# plugins
sevan_plugins := $(patsubst %.C, %-$(METHOD).plugin, $(sevan_plugfiles))
sevan_plugins += $(patsubst %.c, %-$(METHOD).plugin, $(sevan_cplugfiles))
sevan_plugins += $(patsubst %.f, %-$(METHOD).plugin, $(sevan_fplugfiles))
sevan_plugins += $(patsubst %.f90, %-$(METHOD).plugin, $(sevan_f90plugfiles))

# sevan main
sevan_main_src    := $(SEVAN_DIR)/src/main.C
sevan_app_objects := $(patsubst %.C, %.$(obj-suffix), $(sevan_main_src))

# dependency files
sevan_deps := $(patsubst %.C, %.$(obj-suffix).d, $(sevan_srcfiles)) \
              $(patsubst %.c, %.$(obj-suffix).d, $(sevan_csrcfiles)) \
              $(patsubst %.C, %.$(obj-suffix).d, $(sevan_main_src))

# If building shared libs, make the plugins a dependency, otherwise don't.
ifeq ($(libmesh_shared),yes)
  sevan_plugin_deps := $(sevan_plugins)
else
  sevan_plugin_deps :=
endif

all:: $(sevan_LIB)

$(sevan_LIB): $(sevan_objects) $(sevan_plugin_deps)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
	  $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(sevan_objects) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS) -rpath $(SEVAN_DIR)
	@$(libmesh_LIBTOOL) --mode=install --quiet install -c $(sevan_LIB) $(SEVAN_DIR)

# include SEVAN dep files
-include $(sevan_deps)

# how to build SEVAN application
ifeq ($(APPLICATION_NAME),sevan)
all:: sevan

sevan: $(sevan_APP)

$(sevan_APP): $(moose_LIB) $(elk_MODULES) $(sevan_LIB) $(sevan_app_objects)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
          $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(sevan_app_objects) $(sevan_LIB) $(elk_MODULES) $(moose_LIB) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(ADDITIONAL_LIBS)

endif

#
# Maintenance
#
delete_list := $(sevan_APP) $(sevan_LIB) $(SEVAN_DIR)/libsevan-$(METHOD).*

cleanall:: 
	make -C $(SEVAN_DIR) clean 

###############################################################################
# Additional special case targets should be added here
