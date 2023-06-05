
JAR = jar
JAVAC = javac
RM = rm -f

name = SignificanceNoncoding

sources := $(wildcard *.java)

libs = Jama-1.0.3.jar commons-math3-3.6.1.jar
classpath := $(shell echo . $(libs) | tr ' ' ':')

all: $(name).jar

cleanobj:
	$(RM) *.class
	$(RM) manifest.txt

clean: cleanobj
	$(RM) $(name).jar

manifest.txt:
	( echo "Main-Class: $(name)" && \
	  echo "Class-Path: $(libs)" ) >$@

%.class: %.java $(sources) $(libs)
	$(JAVAC) -cp $(classpath) $<

$(name).jar: $(name).class manifest.txt
	$(JAR) -cfm $(name).jar manifest.txt *.class

.PHONY: all cleanobj clean
.SUFFIXES:
