
JAR = jar
JAVAC = javac
CURL = curl
RM = rm -f

name = SignificanceNoncoding

goog_url = https://storage.googleapis.com/noncoding_analysishg19

sources := $(wildcard *.java)

libs = Jama-1.0.3.jar commons-math3-3.6.1.jar
classpath := $(shell echo . $(libs) | tr ' ' ':')

all: $(name).jar

cleanobj:
	$(RM) *.class
	$(RM) manifest.txt

clean: cleanobj
	$(RM) $(name).jar

AnnotationFilesComplete.zip:
	$(CURL) -L -o $@ $(goog_url)/$(@F)

manifest.txt:
	( echo "Main-Class: $(name)" && \
	  echo "Class-Path: $(libs)" ) >$@

%.class: %.java $(sources) $(libs)
	$(JAVAC) -cp $(classpath) $<

$(name).jar: $(name).class manifest.txt
	$(JAR) -cfm $(name).jar manifest.txt *.class

.PHONY: all cleanobj clean
.SUFFIXES:
