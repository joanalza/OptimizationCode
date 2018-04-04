CC = g++ -O2 -Wno-deprecated

all: RankingEDAsCEC

RankingEDAsCEC: main.o PBP.o PFSP.o LOP.o QAP.o TSP.o API.o Tools.o Cayley.o Distance.o Kendall.o GeneralizedKendall.o GeneralizedDistance.o RankingModel.o GeneralizedMallowsModel.o MallowsModel.o Population.o Individual.o RankingEDA.o GeneralizedCayley.o Ferrers_diagram.o Ferrers_diagram2.o Generic.o Ulam.o Ulam2.o NewtonRaphson.o
	$(CC)  main.o PBP.o PFSP.o LOP.o QAP.o TSP.o API.o Tools.o Cayley.o Distance.o Kendall.o GeneralizedKendall.o GeneralizedDistance.o RankingModel.o GeneralizedMallowsModel.o MallowsModel.o Population.o Individual.o RankingEDA.o GeneralizedCayley.o Ferrers_diagram.o Ferrers_diagram2.o Generic.o Ulam.o Ulam2.o NewtonRaphson.o -lstdc++ -lm -o RankingEDAsCEC

main.o: main.cpp RankingEDA.h PBP.h PFSP.h LOP.h QAP.h TSP.h API.h
	$(CC) -c main.cpp
 
Cayley.o: Cayley.cpp Cayley.h Tools.h Distance.h Variables.h NewtonRaphson.h
	$(CC) -c Cayley.cpp

Kendall.o: Kendall.cpp Kendall.h Tools.h Distance.h Variables.h NewtonRaphson.h
	$(CC) -c Kendall.cpp
	
Ulam.o: Ulam.cpp Ulam.h Tools.h Distance.h Variables.h Ferrers_diagram.h NewtonRaphson.h
	$(CC) -c Ulam.cpp
	
Ulam2.o: Ulam2.cpp Ulam2.h Tools.h Distance.h Variables.h Ferrers_diagram2.h NewtonRaphson.h
	$(CC) -c Ulam2.cpp
	
Ferrers_diagram.o: Ferrers_diagram.cpp Ferrers_diagram.h Generic.h
	$(CC) -c Ferrers_diagram.cpp

Ferrers_diagram2.o: Ferrers_diagram2.cpp Ferrers_diagram2.h Generic.h
	$(CC) -c Ferrers_diagram2.cpp
	
Generic.o: Generic.cpp Generic.h
	$(CC) -c Generic.cpp

Distance.o: Distance.cpp Distance.h 
	$(CC) -c Distance.cpp

GeneralizedKendall.o: GeneralizedKendall.cpp GeneralizedKendall.h Tools.h GeneralizedDistance.h NewtonRaphson.h
	$(CC) -c GeneralizedKendall.cpp

GeneralizedCayley.o: GeneralizedCayley.cpp GeneralizedCayley.h Tools.h GeneralizedDistance.h NewtonRaphson.h
	$(CC) -c GeneralizedCayley.cpp

GeneralizedDistance.o: GeneralizedDistance.cpp GeneralizedDistance.h
	$(CC) -c GeneralizedDistance.cpp

RankingModel.o: RankingModel.cpp RankingModel.h Population.h
	$(CC) -c RankingModel.cpp

MallowsModel.o: MallowsModel.cpp MallowsModel.h Individual.h RankingModel.h Population.h Distance.h Kendall.h Cayley.h Tools.h Ulam.h Ulam2.h
	$(CC) -c MallowsModel.cpp

GeneralizedMallowsModel.o: GeneralizedMallowsModel.cpp GeneralizedMallowsModel.h Individual.h RankingModel.h Population.h GeneralizedCayley.h  GeneralizedDistance.h GeneralizedKendall.h Tools.h
	$(CC) -c GeneralizedMallowsModel.cpp

Population.o: Population.cc Population.h Individual.h
	$(CC) -c Population.cc

Individual.o: Individual.cc Individual.h Variables.h
	$(CC) -c Individual.cc

RankingEDA.o: RankingEDA.cpp RankingEDA.h MallowsModel.h GeneralizedMallowsModel.h Tools.h PBP.h Population.h RankingModel.h
	$(CC) -c RankingEDA.cpp

Tools.o: Tools.cpp Tools.h
	$(CC) -c Tools.cpp

PBP.o: PBP.cpp PBP.h
	$(CC) -c PBP.cpp

PFSP.o: PFSP.cpp PFSP.h Tools.h PBP.h
	$(CC) -c PFSP.cpp

LOP.o: LOP.cpp LOP.h Tools.h PBP.h                        
	$(CC) -c LOP.cpp

QAP.o: QAP.cpp QAP.h Tools.h PBP.h                        
	$(CC) -c QAP.cpp

TSP.o: TSP.cpp TSP.h Tools.h PBP.h                        
	$(CC) -c TSP.cpp

API.o: API.cpp API.h Tools.h PBP.h                        
	$(CC) -c API.cpp
	
NewtonRaphson.o: NewtonRaphson.cpp NewtonRaphson.h                         
	$(CC) -c NewtonRaphson.cpp

clean:
	rm -f core *.o RankingEDAsCEC

