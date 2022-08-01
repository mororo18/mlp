CC = gcc
BIN_NAME = solve

MODE = MATRIX

COPT = -O3 -std=c99 -D$(MODE) $(Wno) $(W) -g
W = -Winline
Wno = -Wno-unused-result 
LOPT = -lm

OBJDIR = obj

SRCS = $(wildcard *.c)
OBJS = $(patsubst %.c, $(OBJDIR)/%.o, $(SRCS))



$(BIN_NAME): $(OBJS)
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CC) $^ $(LOPT) -o $(BIN_NAME)


$(OBJDIR)/%.o: %.c
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CC) -c $^ -o $@ $(COPT) 
	
	
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CC) -MM $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

flat: MODE = FLAT
flat: $(BIN_NAME)

clean:
	@echo "\033[31mcleaning obj directory \033[0m"
	@rm $(OBJDIR)/* $(BIN_NAME)

rebuild: clean $(BIN_NAME)

run: $(BIN_NAME)
	@./$(BIN_NAME)
