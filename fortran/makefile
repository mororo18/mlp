COMP = gfortran
BIN_NAME = solve

PREPROC = -cpp 
PREPROC += -DDOUBLE
COPT =  $(W) -g -Ofast
W = -Winline

OBJDIR = obj
MODDIR = mod

SRCS = $(wildcard *.f90)
OBJS = $(patsubst %.f90, $(OBJDIR)/%.o, $(SRCS))

$(BIN_NAME): $(OBJS)
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(COMP) $^ -o $(BIN_NAME) -g
	@rm -f *.mod ||:


$(OBJDIR)/%.o: %.f90
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(COMP) $(PREPROC) $(COPT) -c $< -o $@ -I$(MODDIR) -J$(MODDIR)
	
	
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(COMP) $(PREPROC) $(COPT) -I$(MODDIR) -MM $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

# dependencias
$(OBJDIR)/main.o: $(filter-out $(OBJDIR)/main.o,$(OBJS))
$(OBJDIR)/data.o: $(OBJDIR)/types.o

debug: PREPROC += -DDEBUG 
debug: COPT += -fbounds-check
debug: $(BIN_NAME)

clean:
	@echo "\033[31mcleaning obj directory \033[0m"
	@rm -f $(OBJDIR)/* $(BIN_NAME) $(MODDIR)/* *.o *.mod

rebuild: clean $(BIN_NAME)

run: $(BIN_NAME)
	@./$(BIN_NAME)
