FC       = gfortran
FFLAGS   = -O2 -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow
LDFLAGS  =

MOD_DIR  = -J ../../include/test
BIN_DIR  = ../../bin
SRC_DIR  = ../../src
OBJ_DIR  = ../../obj/test
BIN_LIST = test
TARGET   = $(addprefix $(BIN_DIR)/, $(BIN_LIST))
SRC_LIST = util.f90 debug.f90 io.f90 element.f90 matrix.f90 update.f90 solver.f90 analysis.f90 test.f90
SOURCES  = $(addprefix $(SRC_DIR)/, $(SRC_LIST))
OBJS     = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))
RM       = rm

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(MOD_DIR) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET) ../../include/*.mod

.PHONY: clean
