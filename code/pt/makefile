# Compiler settings
CXX=g++
CXXFLAGS=-std=c++11 -Wall

# Directories
SRCDIR=src
OBJDIR=obj
BINDIR=bin

# Target executable name
TARGET=$(BINDIR)/program

# Source and object files
SOURCES=$(wildcard $(SRCDIR)/*.cpp)
OBJECTS=$(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# Default target
default: $(TARGET)

# Rule to create the target directory
$(BINDIR):
	mkdir -p $(BINDIR)

# Rule to create the object directory
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Linking the target binary
$(TARGET): $(OBJECTS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

# Compiling source files to object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -rf $(OBJDIR) $(BINDIR)
