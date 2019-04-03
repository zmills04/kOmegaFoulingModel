#if !defined(AFX_REDUCETESTS_H__INCLUDED_)
#define AFX_REDUCETESTS_H__INCLUDED_

#include "ReduceGenerator.h"
#include "Reducer.h"



namespace Tests
{
	template <typename T>
	void printNameAndValue(ReducerBase<T> &rw_, double actualval)
	{
		std::cout << "Reducer Name: " << rw_.Name << ", Value: " << std::setprecision(7) <<
			rw_.reduceSingle() << ", Actual Value: " << std::setprecision(7) << actualval << std::endl;
	}

	// TODO: test that the kernel is not accessing an out of bounds element by
	// creating a kernel and reducing a subsection of that kernel
	void testReduce()
	{
		int FullSize = 1048576;
		int FullSizeX = 4096, FullSizeY = 256;
		Array1Dd A1, A2, A3, A4, A5, A6;
		Array2Dd B1, B2, B3, B4;

		A1.zeros(FullSize);
		A1.setName("A1");
		A1.fill(1.);
		A1.allocate_buffer_w_copy();
		double A1val = 1048576.;



		A2.setName("A2");
		A2.zeros(FullSize - 54, FullSize);
		A2.fill(1.);
		A2.allocate_buffer_w_copy();
		double A2val = 1048576 - 54.;


		A3.setName("A3");
		A3.zeros(FullSize - WORKGROUPSIZE_RED + 10);
		A3.fill(1.);
		A3.allocate_buffer_size(FullSize);
		A3.copy_to_buffer_size(FullSize);
		double A3val = 1048576. - (double)WORKGROUPSIZE_RED + 10.;

		A4.setName("A4");
		A4.zeros(FullSize - 108);
		A4.fill(1.);
		A4.allocate_buffer_w_copy();
		double A4val = 1048576. - 108.;

		B1.setName("B1");
		B1.zeros(FullSizeX, FullSizeY);
		B1.fill(1.);
		B1.allocate_buffer_w_copy();
		double B1val = (double)(FullSizeX*FullSizeY);


		B2.setName("B2");
		B2.zeros(FullSizeX - 54, FullSizeX, FullSizeY - 12, FullSizeY);
		B2.fill(1.);
		B2.allocate_buffer_w_copy();
		double B2val = (double)((FullSizeX - 54)*(FullSizeY - 12));

		B3.setName("B3");
		B3.zeros(FullSizeX - WORKGROUPSIZE_RED + 10, FullSizeY - 10);
		B3.fill(1.);
		B3.allocate_buffer_size(FullSizeX*FullSizeY);
		B3.copy_to_buffer();
		double B3val = (double)((FullSizeX - WORKGROUPSIZE_RED + 10)*(FullSizeY - 10));


		B4.setName("B4");
		B4.zeros(FullSizeX - 108, FullSizeY - 33);
		B4.fill(1.);
		B4.allocate_buffer_w_copy();
		double B4val = (double)((FullSizeX - 108)*(FullSizeY - 33));

		A5.setName("A5");
		A5.zeros(FullSize + 4343);
		A5.fill(1.);
		A5.allocate_buffer_w_copy();
		double A5val = 1048576. + 4343.;

		A6.setName("A6");
		A6.zeros(2 * FullSize);
		A6.fill(1.);
		A6.allocate_buffer_w_copy();
		A6.resetSizesDebug(FullSize);
		double A6val = A1val;



		Reducer<double, ReduceGenerator::Sum, 100> redA1, redA2,
			redA3, redA4, redA5, redA6, redB1, redB2, redB3, redB4;
		Reducer<double, ReduceGenerator::Sqr, 100> redSqrA1;
		Reducer<double, ReduceGenerator::Dot, 100> redDotA1A2;
		Reducer<double, ReduceGenerator::Abs, 100> redAbsA1;
		// Need to figure out why work group size array is remaining null after being filled
		redDotA1A2.ini(A1, A2, "reduceDotA1A2"); // working
		redA1.ini(A1, "reduceA1"); // working
		redA2.ini(A2, "reduceA2"); // working
		redA3.ini(A3, "reduceA3"); // working
		redA4.ini(A4, "reduceA4"); // working
		redA5.ini(A5, "reduceA5"); // working
		redA6.ini(A6, "reduceA6"); // working
		redB1.ini(B1, "reduceB1"); // working
		redB2.ini(B2, "reduceB2"); // working
		redB3.ini(B3, "reduceB3"); // working
		redB4.ini(B4, "reduceB4"); // working

		redSqrA1.ini(A1, "reduceSqrA1"); // working
		redAbsA1.ini(A1, "reduceAbsA1"); // working

		ReduceGenerator::ReduceInstance()->buildSource();

		printNameAndValue(redA1, A1val);
		printNameAndValue(redA2, A2val);
		printNameAndValue(redA3, A3val);
		printNameAndValue(redA4, A4val);
		printNameAndValue(redA5, A5val);
		printNameAndValue(redA6, A6val);
		printNameAndValue(redB1, B1val);
		printNameAndValue(redB2, B2val);
		printNameAndValue(redB3, B3val);
		printNameAndValue(redB4, B4val);

		A1.fill(2.);
		A1.copy_to_buffer();
		printNameAndValue(redSqrA1, A1val * 4);

		A2.fill(3.);
		A2.copy_to_buffer();
		printNameAndValue(redDotA1A2, 6.*MIN(A1val, A2val));

		A1.fill(-1.);
		A1.copy_to_buffer();
		printNameAndValue(redAbsA1, A1val);

	}
}
#endif