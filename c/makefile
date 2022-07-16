#detecta se o sistema é de 32 ou 64 bits
BITS_OPTION = -m64

#### define o compilador
#CPPC = clang++
CPPC = gcc
#############################

#### opcoes de compilacao e includes
CCOPT = $(BITS_OPTION) -fPIC -O3 -fexceptions -DIL_STD -std=c99 -lm -g
CONCERTINCDIR = $(CONCERTDIR)/include
#############################

#### flags do linker
#############################

#### diretorios com os source files e com os objs files
#SRCDIR = src
OBJDIR = obj
#############################

#### lista de todos os srcs e todos os objs
SRCS = $(wildcard *.c)
#SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst %.c, $(OBJDIR)/%.o, $(SRCS))
#OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))
#############################

#### regra principal, gera o executavel
solve: $(OBJS) 
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) -o $@ $(CPPFLAGS) $(CCOPT) -g 
############################

#inclui os arquivos de dependencias
-include $(OBJS:.o=.d)

#regra para cada arquivo objeto: compila e gera o arquivo de dependencias do arquivo objeto
#cada arquivo objeto depende do .c e dos headers (informacao dos header esta no arquivo de dependencias gerado pelo compiler)
$(OBJDIR)/%.o: %.c
#$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CPPFLAGS) $(CCOPT) -c $< -o $@ -g
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC)   -MM $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

#delete objetos e arquivos de dependencia
clean:
	@echo "\033[31mcleaning obj directory \033[0m"
	@rm tsp -f $(OBJDIR)/*.o $(OBJDIR)/*.d


rebuild: clean tsp
