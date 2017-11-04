exec    = MC.exe
cc      = g++
opt     = -O3
cflags  = -std=c++11 $(opt)

main    = Main.cpp
objects = $(patsubst %.cpp,%.o,$(filter-out $(main), $(wildcard *.cpp)))

.PHONY : all clean

all :	$(objects) 
	@rm -f $(exec)
	@$(MAKE) $(exec)

%.o : %.cpp
	$(cc) $(cflags) -c $<

$(exec) : $(main)
	$(cc) $(cflags) $(objects) $< -o $@

clean :
	rm -f $(objects) $(exec)
