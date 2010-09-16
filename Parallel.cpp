#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <dps/dps.h>

#define MAX_IT 15
#define EPSI 0.01

/***********************************
 *
 *		Class Vector
 *
 ***********************************/

class Vector : public dps::AutoSerial
{
	CLASSDEF(Vector)
	MEMBERS
		ITEM(int,size)
		ITEM(dps::Buffer<double>,data)
	CLASSEND;

public:
	Vector() {}

	Vector(int s)
	{
		this->size = s;
		data.resize(size);
	}

	Vector(int s, double value)
	{
		this->size = s;
		data.resize(size);
		for(int i=0; i<s; i++)
		{
			data[i] = value;
		}
	}

	Vector(dps::Buffer<double> array1D)
	{
		this->size = (int)array1D.size();
		data = array1D;
	}

	Vector(double* array1D, int s)
	{
		this->size = s;
		data = dps::Buffer<double>(array1D, s);
	}

	inline double& operator()(int i)
	{
		return data[i];
	}

	Vector operator-(Vector& V)
	{
		Vector R = Vector(size);
		for(int i=0; i<V.size; i++)
		{
			R.data[i] = data[i] - V.data[i];
		}
		return R;
	}

	std::string toString()
	{
		std::ostringstream o;
		for(int i=0; i<size; i++)
		{
			o << data[i] << "\t";
		}
		return o.str();
	}

	double sumabs()
	{
		double sabs = 0;
		for(int i=0; i<size; i++)
		{
			sabs += fabs(data[i]);
		}
		return sabs;
	}

	static Vector generateRandomB(int s)
	{
		Vector V = Vector(s);
		srand( 3687687 );
		for(int i=0; i<V.size; i++)
		{
			V.data[i] = rand() - RAND_MAX/2;
		}
		return V;
	}

	static Vector zeros(int s)
	{
		Vector V = Vector(s);
		for(int i=0; i<V.size; i++)
		{
			V.data[i] = 0;
		}
		return V;
	}

	Vector get(int firstIndex, int lineCount)
	{
		return Vector(&data[firstIndex], lineCount);
	}

};



/***********************************
 *
 *		Class Matrix
 *
 ***********************************/

class Matrix : public dps::AutoSerial
{
	CLASSDEF(Matrix)
	MEMBERS
		ITEM(int,width)
		ITEM(int,height)
		ITEM(dps::Buffer<double>,data)
	CLASSEND;

public:
	Matrix() {}

	Matrix(int h, int w)
	{
		this->width = w;
		this->height = h;
		data.resize(w*h);
	}

	Matrix(dps::Buffer<double> array1D, int w)
	{
		this->width = w;
		this->height = (int)array1D.size()/w;
		data = array1D;
	}

	Matrix(double* array1D, int h, int w)
	{
		this->width = w;
		this->height = h;
		data = dps::Buffer<double>(array1D, w*h);
	}

	inline double& operator()(int i, int j)
	{
		return data[i*width+j];
	}

	inline Vector operator()(int i)
	{
		return Vector(&data[i],width);
	}

	std::string toString()
	{
		std::ostringstream o;
		for(int i=0; i<height; i++)
		{
			int ov = i*width;
			for(int j=0; j<width; j++)
			{
				o << data[j+ov] << "\t";
			}
			o << "\n";
		}
		return o.str();
	}

	static Matrix generateRandomA(int h, int w)
	{
		Matrix M = Matrix(h,w);
		srand(2834863);
		for(int i=0; i<h; i++)
		{
			double lineSum = 0;
			int ov = i*w;
			for(int j=0; j<w; j++)
			{
				if(i!=j)
				{
					// generate a number between -RAND_MAX/2 and RAND_MAX/2;
					double v = rand() - RAND_MAX/2;
					M.data[ov+j] = v;
					lineSum += abs(v);
				}
			}
			double r = rand();
			int posOrNeg = (r-RAND_MAX/2 >= 0)?1:-1;
			M.data[ov+i] = posOrNeg * ( lineSum + r );
		}
		return M;
	}

	static Matrix zeros(int h, int w)
	{
		Matrix M = Matrix(h,w);
		for(int i=0; i<h; i++)
		{
			int ov = i*w;
			for(int j=0; j<w; j++)
			{
				M.data[ov+j] = 0;
			}
		}
		return M;
	}

	Matrix getLines(int firstLine, int lineCount)
	{
		return Matrix(&data[firstLine*width], lineCount, width);
	}

};
/*
Vector operator*(Matrix& A, Vector& B)
{
	Vector R = Vector(height);
	for(int i=0; i<A.height; i++)
	{
		double lineprod = 0;
		for(int j=0; j<A.width; j++)
		{
			lineprod += A(i,j) * B(j);
		}
		R(i) = lineprod;
	}
	return R;
}
*/








/***********************************
 ***********************************
 *
 *     Paralle Implementation
 *
 ***********************************
 ***********************************/




class SplitInInitialDataObject : public dps::AutoSerial
{
	CLASSDEF(SplitInInitialDataObject)
	MEMBERS
		ITEM(Matrix,A)
		ITEM(Vector,b)
		ITEM(Int32,X0)
	CLASSEND;
};

class SplitOutInitialDataObject : public dps::AutoSerial
{
	CLASSDEF(SplitOutInitialDataObject)
	MEMBERS
		ITEM(Matrix,Apart)
		ITEM(Vector,bpart)
		ITEM(Int32,X0)
		ITEM(Int32,firstLine)
		ITEM(Int32,lineCount)
		ITEM(Int32,Aheight)
	CLASSEND;

	SplitOutInitialDataObject(Int32 firstLine = 0, Int32 lineCount = 0) 
	{
		this->firstLine = firstLine;
		this->lineCount = lineCount; 
	}

	SplitOutInitialDataObject(SplitInInitialDataObject* in, Int32 firstLine, Int32 lineCount) 
	{
		this->Apart = in->A.getLines(firstLine, lineCount);
		this->bpart = in->b.get(firstLine, lineCount);
		this->X0 = in->X0;
		this->Aheight = in->A.height;
		this->firstLine = firstLine;
		this->lineCount = lineCount; 
	}
};


class MergeInitialDataObject : public dps::AutoSerial
{
	CLASSDEF(MergeInitialDataObject)
	MEMBERS
		ITEM(bool,ok)
	CLASSEND;
};

class RoundDataObject : public dps::AutoSerial
{
	CLASSDEF(RoundDataObject)
	MEMBERS
		ITEM(Vector,Xk)
		ITEM(Int32,k)
		ITEM(Int32,firstLine)
		ITEM(Int32,lineCount)
	CLASSEND;

	RoundDataObject(Int32 firstLine = 0, Int32 lineCount = 0) 
	{
		this->firstLine = firstLine;
		this->lineCount = lineCount; 
	}
};

class ProcessThread
{
	IDENTIFY(ProcessThread);
public:
	Matrix Apart;
	Vector bpart;
	Int32 X0;
	Int32 firstLine;
	Int32 lineCount;
	Int32 Aheight;
	virtual ~ProcessThread() { }

};


class MainThread
{
	IDENTIFY(MainThread);
public:
	Int32 slaveCount;
	Int32 Aheight;
	Double startTime;
	Int32 maxIt;
};

class SplitInitial : public dps::SplitOperation<SplitInInitialDataObject,SplitOutInitialDataObject,MainThread>
{
	IDENTIFY(SplitInitial);

	void execute(SplitInInitialDataObject* in)
	{
		//dps::DPSLog.write(0) << "SplitInitial: sending initial data to nodes...";
		getThread()->slaveCount = getController()->getThreadCollection<ProcessThread>("process").getSize();
		//dps::DPSLog.write(0) << "SplitInitial: number of slaves : " << getThread()->slaveCount;
		getThread()->Aheight = in->A.height;
		int nodeCount = getThread()->slaveCount;
		int linePerNode = in->A.height/nodeCount;
		int lineRest = in->A.height % nodeCount;
		//dps::DPSLog.write(0) << "SplitInitial: line per slave : " << linePerNode << ". Extra amount of lines : " << lineRest;
		int startLine = 0;
		// give the extra amount of lines to the first node
		int endLine = linePerNode + lineRest;
		int i = 1;
		while(startLine < in->A.height)
		{
			//dps::DPSLog.write(0) << "SplitInitial: slave " << i << " startline:" << startLine << " linecount:" << endLine-startLine;
			SplitOutInitialDataObject* out = new SplitOutInitialDataObject(in, startLine, endLine-startLine);
			postDataObject(out);
			startLine = endLine;
			endLine += linePerNode;
			//dps::DPSLog.write(0) << "SplitInitial: data sent to slave " << i++;
		}
		//dps::DPSLog.write(0) << " SplitInitial: terminated!" << std::endl;
	}
	
};

class StoreInitial : public dps::LeafOperation<SplitOutInitialDataObject,MergeInitialDataObject,ProcessThread>
{
	IDENTIFY(StoreInitial);

public:
	void execute(SplitOutInitialDataObject *in)
	{
		//dps::DPSLog.write(0) << "StoreInitial " << in->firstLine/in->lineCount;
		getThread()->Apart = in->Apart;
		getThread()->bpart = in->bpart;
		getThread()->X0 = in->X0;
		getThread()->firstLine = in->firstLine;
		getThread()->lineCount = in->lineCount;
		getThread()->Aheight = in->Aheight;
		postDataObject(new MergeInitialDataObject());
	}
};

class MergeInitial : public dps::MergeOperation<MergeInitialDataObject,RoundDataObject,MainThread>
{
	IDENTIFY(MergeInitial);

	void execute(MergeInitialDataObject* in)
	{
		while(waitForNextDataObject()!=NULL);
		//dps::DPSLog.write(0) << "MergeInitial: starting...";
		RoundDataObject* firstRound = new RoundDataObject();
		firstRound->k = 0;
		postDataObject(firstRound);
		//dps::DPSLog.write(0) << "MergeInitial: terminated!" << std::endl;
	}
};


class DistributeRound : public dps::SplitOperation<RoundDataObject,RoundDataObject,MainThread>
{
	IDENTIFY(DistributeRound);

	void execute(RoundDataObject* in)
	{
		//dps::DPSLog.write(0) << "Round " << in->k;
		//dps::DPSLog.write(0) << "DitributeRound: starting...";
		int nodeCount = getThread()->slaveCount;
		//dps::DPSLog.write(0) << "DitributeRound: " << nodeCount << " slaves!";
		for(int i=0; i<nodeCount; i++)
		{
			int firstLine = i;
			int lineCount = 1;
			RoundDataObject* out = new RoundDataObject(firstLine, lineCount);
			out->k = in->k + 1;
			if(out->k > 1) out->Xk = in->Xk;
			postDataObject(out);
		}
		//dps::DPSLog.write(0) << "DitributeRound: terminated!";
	}
};

class ProcessRound : public dps::LeafOperation<RoundDataObject,RoundDataObject,ProcessThread>
{
	IDENTIFY(ProcessRound);
	
public:
	static Vector jacobiRound(Matrix& Apart, Vector& bpart, Vector& Xk, int firstLine, int lineCount)
	{
		Vector Xk1part = Vector(lineCount);
		int columns = Apart.width;
		for(int i=0; i<lineCount; i++)
		{
			double s = 0;
			int ipart = i+firstLine;
			for(int j=0; j<ipart; j++)
			{
				s += Apart(i,j) * Xk(j);
			}
			for(int j=ipart+1; j<columns; j++)
			{
				s += Apart(i,j) * Xk(j);
			}
			Xk1part(i) = ( bpart(i) - s ) / Apart(i,ipart);
		}
		return Xk1part;
	}

	void execute(RoundDataObject *in)
	{
		//dps::DPSLog.write(0) << "ProcessRound" << in->firstLine/in->lineCount << ": starting...";
		Int32 firstLine = getThread()->firstLine;
		Int32 lineCount = getThread()->lineCount;
		if(in->k == 1) in->Xk = Vector(getThread()->Aheight, getThread()->X0);
		RoundDataObject* out = new RoundDataObject(firstLine, lineCount);
		out->Xk = jacobiRound(getThread()->Apart, getThread()->bpart, in->Xk, firstLine, lineCount);
		out->firstLine = firstLine;
		out->lineCount = lineCount;
		out->k = in->k;
		postDataObject(out);
	}
};

class MergeRound : public dps::MergeOperation<RoundDataObject,RoundDataObject,MainThread>
{
	IDENTIFY(MergeRound);

	void execute(RoundDataObject *in)
	{
		//dps::DPSLog.write(0) << "MergeRound: starting...";
		RoundDataObject *rdo = new RoundDataObject();
		rdo->Xk = Vector(getThread()->Aheight);
		rdo->k = in->k;
		do
		{
			int i = 0;
			//dps::DPSLog.write(0) << "MergeRound: in->Xk=" << in->Xk.toString();
			while(i < in->lineCount)
			{
				int ipart = i + in->firstLine;
				rdo->Xk(ipart) = in->Xk(i);
				i++;
			}
		} while((in=waitForNextDataObject())!=NULL);
		//dps::DPSLog.write(0) << "MergeRound: Xk=" << rdo->Xk.toString();
		postDataObject(rdo);
		//dps::DPSLog.write(0) << "MergeRound: terminated!";
	}
};


class HasConverged : public dps::Condition<RoundDataObject>
{
	IDENTIFY(HasConverged);
	Vector Xkprev;

public:
	bool condition(RoundDataObject *in)
	{
		bool loop = true;
		if(in->k >= MAX_IT)
		{
			return !loop;
		}
		else if(in->k > 1)
		{
			double ununm1 = (in->Xk - Xkprev).sumabs();
			//dps::DPSLog.write(0) << "HasConverged: ununm1=" << ununm1;
			loop = (ununm1 > EPSI);
		}
		Xkprev = in->Xk;
		return loop;
	}
};

//! Round robin routing function
class MyRoundRobinRoute : public dps::Route<RoundDataObject>
{
	IDENTIFY(MyRoundRobinRoute);
public:
	//! Route according to the part of the matrix
	virtual Size route(RoundDataObject *in) { return in->firstLine/in->lineCount; }
};

//! Round robin routing function
class MyInitialRoute : public dps::Route<SplitOutInitialDataObject>
{
	IDENTIFY(MyInitialRoute);
public:
	//! Route according to the part of the matrix
	virtual Size route(SplitOutInitialDataObject *in){ return in->firstLine/in->lineCount; }
};

//! Constant routing function
class MainRoute : public dps::Route<dps::ISerializable>
{
	IDENTIFY(MainRoute);
public:
	//! Route all data objects to thread 0
	virtual Size route(dps::ISerializable *in) { return 0; }
};


//! Application class
class JacobiApp : public dps::Application
{
public:
	virtual bool init();
	virtual void start();
	Vector serial(Matrix& A, Vector& b, double x0);
	void bench(dps::Flowgraph& theGraph);

	IDENTIFY(JacobiApp);
};

bool JacobiApp::init()
{
	return true;
}

Vector JacobiApp::serial(Matrix& A, Vector& b, double x0)
{
	int lines = A.height;
	int columns = A.width;
	Vector Xk = Vector(lines,x0);
	Vector Xk1 = Vector(lines);
    int k = 0;
	while(k < MAX_IT)
	{
		for(int i=0; i<lines; i++)
		{
			double s = 0;
			for(int j=0; j<i; j++)
			{
				s += A(i,j) * Xk(j);
			}
			for(int j=i+1; j<columns; j++)
			{
				s += A(i,j) * Xk(j);
			}
			Xk1(i) = ( b(i) - s ) / A(i,i);
		}
		//dps::DPSLog.write(0) << "k = " << k;
		//dps::DPSLog.write(0) << "b = " << b.toString() << " \ns = " << s1.toString();
		double ununm1 = (Xk1-Xk).sumabs();
		//double aunmb = ((A*Xk1)-b).sumabs();
		//printf("it: %d conv1: %f conv2: %f\n", k, ununm1, aunmb);
		Xk = Xk1;
		if(ununm1<=EPSI) break;
		k++;
	}
	return Xk;
}

void JacobiApp::bench(dps::Flowgraph& theGraph)
{
	double st,et;
	std::ofstream myfile;
	myfile.open("results.csv");
	myfile << "size;parallel;serial\n";
	for(int rsize=10; rsize<=10000; (rsize>=100)?rsize+=100:rsize+=10)
	{
		SplitInInitialDataObject *initial = new SplitInInitialDataObject();
		initial->b = Vector::generateRandomB(rsize);
		initial->X0 = 1;
		initial->A = Matrix::generateRandomA(rsize,rsize);

		//parallel
		st=dps::Timer::getMillis();
		RoundDataObject *result = (RoundDataObject*)getController()->callSchedule(theGraph,initial);
		et=dps::Timer::getMillis();
		double parallelTime = et-st;

		//serial
		st=dps::Timer::getMillis();
		Vector X = serial(initial->A, initial->b, initial->X0);
		et=dps::Timer::getMillis();
		double serialTime = et-st;
		myfile << rsize << ";" << parallelTime << ";" << serialTime << "\n";
		delete initial;
		result->release();
	}
	myfile.close();
}

void JacobiApp::start()
{
	dps::ThreadCollection<MainThread> theMainThread = 
		getController()->createThreadCollection<MainThread>("main");
	dps::ThreadCollection<ProcessThread> processThreads = 
		getController()->createThreadCollection<ProcessThread>("process",1); //2);

	if(DPS_FAILED(theMainThread.addThread("0")))
	{
		std::cout << "Could not map main thread collection" << std::endl;
		return;
	}
	const char *pattern=getController()->getConfig().getValue("pat","2");
	if(DPS_FAILED(processThreads.addThread(dps::MPIMapper::get(pattern))))
	{
		std::cout << "Could not map process thread collection" << std::endl;
		return;
	}

	dps::FlowgraphNode<DistributeRound,MainRoute> split(theMainThread);
	
	dps::FlowgraphBuilder theGraphB =
		dps::FlowgraphNode<SplitInitial,MainRoute>(theMainThread) >>
		dps::FlowgraphNode<StoreInitial,MyInitialRoute>(processThreads) >>
		dps::FlowgraphNode<MergeInitial,MainRoute> (theMainThread) >>
		split >>
		dps::FlowgraphNode<ProcessRound,MyRoundRobinRoute>(processThreads) >>
		dps::FlowgraphNode<MergeRound,MainRoute>(theMainThread) >>
		dps::FlowgraphLoop<HasConverged>(split);

	dps::Flowgraph theGraph=getController()->createFlowgraph("graph", theGraphB);
	
	bench(theGraph);

	/*
	int rsize = atoi(getController()->getConfig().getValue("size","5"));
	dps::DPSLog.write(0) << "Size of the Matrix A is " << rsize << " x " << rsize;
	
	SplitInInitialDataObject *initial = new SplitInInitialDataObject();
	initial->b = Vector::generateRandomB(rsize);
	initial->X0 = 1;
	initial->A = Matrix::generateRandomA(rsize,rsize);
	
	dps::DPSLog.write(0) << "A = " << initial->A.toString();
	dps::DPSLog.write(0) << "b = " << initial->b.toString();
	dps::DPSLog.write(0) << "X0(0) = " << initial->X0;
	
	dps::DPSLog.write(0) << "Starting parallel execution";
	double st=dps::Timer::getMillis();
	RoundDataObject *result = (RoundDataObject*)getController()->callSchedule(theGraph,initial);
	double et=dps::Timer::getMillis();
	dps::DPSLog.write(0) << "Done in " << (et-st) << " ms \nX = " << result->Xk.toString();
	result->release();
	*/
}



//! Starts up application
int main(int argc, char *argv[])
{
	return dps::dpsMain(argc,argv,new JacobiApp());
}

