#include "Test.hpp"

#ifdef COMPILE_TEST



void test_OnDiskMatrix() {
	OnDiskMatrix<double> matrix{ "test.mat",10,10 };
	// write rows into the matrix
	for (auto i = 0; i < 10; ++i) {
		auto row = Eigen::MatrixXd(1, 10);
		row.setOnes();
		row *= (double)i / 9.0;
		cout << fixed << setprecision(4) << row << "\n";
		matrix.write_row(row, i);
	}
	cout << "\n";

	for (auto i = 0; i < 10; ++i) {
		cout << fixed << setprecision(4) << matrix.read_row(i) << "\n";
	}
	cout << "\n";

	matrix.generate_transpose_matrix("trans.mat");

	OnDiskMatrix<double> trans{ "trans.mat" };
	for (auto i = 0; i < 10; ++i) {
		cout << fixed << setprecision(4) << trans.read_row(i) << "\n";
	}
}

void test_GenerateRandomMatrix() {
	OnDiskMatrix<double> ondisk{ "random.mat",5,10 };
	for (auto i = 0; i < 5; ++i) {
		ondisk.write_row(Eigen::MatrixXd::Random(1, 10), i);
	}
}

void test_OnDiskMatrix_ReadingTime() {
	Timer timer;
	srand(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
	const int rows = 2000;
	const int cols = 2000;

	auto matrix_size = (double)(rows*cols * sizeof(double)) / (double)(1024 * 1024);

	cout << "constructing matrix...\n";
	timer.begin_timing();

	Eigen::MatrixXd matrix{ rows,cols };
	matrix.setRandom();
	timer.stop_timing();
	cout << "matrix constructed in " << timer.get_duration() << " seconds.\n";
	cout << "the size of matrix is " << matrix_size << " MB.\n\n";

	cout << "saving matrix into file...\n";
	timer.begin_timing();
	OnDiskMatrix<double> disk_matrix{ "test.mat",rows,cols };
	for (auto i = 0; i < rows; ++i) {
		disk_matrix.write_row(matrix.row(i), i);
	}
	timer.stop_timing();
	cout << "matrix saved in " << timer.get_duration() << " seconds.\n\n";

	cout << "test saving transpose matrix into file...\n";
	timer.begin_timing();
	disk_matrix.generate_transpose_matrix("trans.mat");
	timer.stop_timing();
	cout << "finished in " << timer.get_duration() << " seconds.\n\n";

	Eigen::MatrixXd row{ 1,cols };
	cout << "conducting row reading test in RAM...\n";
	timer.begin_timing();
	for (auto i = 0; i < rows; ++i) {
		row = move(matrix.row(i));
	}
	timer.stop_timing();
	cout << "reading finished in " << timer.get_duration() << " seconds.\n\n";

	cout << "conducting row reading test from disk...\n";
	timer.begin_timing();
	for (auto i = 0; i < rows; ++i) {
		row = move(disk_matrix.read_row(i));
	}
	timer.stop_timing();
	cout << "reading finished in " << timer.get_duration() << " seconds.\n\n";

	auto get_random_pos = [&]() {
		auto pos1 = rand() % rows;
		auto pos2 = rand() % cols;
		return pair<int, int>{pos1, pos2};
	};
	cout << "conducting random accessing test...(50000 elements)\n";
	timer.begin_timing();
	for (auto i = 0; i < 50000; ++i) {
		auto pos = get_random_pos();
		auto ele = matrix(pos.first,pos.second);
	}
	timer.stop_timing();
	cout << "RAM matrix finished in " << timer.get_duration() << " seconds.\n";

	timer.begin_timing();
	for (auto i = 0; i < 5000; ++i) {
		auto pos = get_random_pos();
		auto ele = disk_matrix.get_element(pos.first, pos.second);
	}
	timer.stop_timing();
	cout << "disk matrix finished in " << timer.get_duration() << " seconds.\n\n";

	cout << "running row equality test...\n";
	timer.begin_timing();
	for (auto i = 0; i < rows; ++i) {
		auto ram_row = move(matrix.row(i));
		auto disk_row = move(disk_matrix.read_row(i));
		for (auto j = 0; j < cols; ++j) {
			expr_check(abs(ram_row(0, j) - disk_row(0, j)) < 1.0e-10, "elements are different");
		}
	}

	timer.stop_timing();
	cout << "equality test finished in " << timer.get_duration() << "seconds.\n\n";

	cout << "running random access equality test...(10000 elements)\n";
	timer.begin_timing();
	for (auto i = 0; i < 10000; ++i) {
		auto pos = get_random_pos();
		expr_check(abs(matrix(pos.first, pos.second) - disk_matrix.get_element(pos.first, pos.second)) < 1.0e-10, "elements are different");
	}
	timer.stop_timing();
	cout << "finished in " << timer.get_duration() << " seconds.\n\n";

}

void test_SimplexMethod() {
	Eigen::MatrixXd mat{ 3,5 };
	mat << 1., -2., 1., 1., 0.,
		-4., 1., 2., 0., -1.,
		-2., 0., 1., 0., 0.;
	Eigen::MatrixXd vec_b{ 3,1 };
	vec_b << 11., 3., 1.;
	Eigen::MatrixXd vec_c{ 1,5 };
	vec_c << 3., -1., -1., 0., 0.;

	OnDiskMatrix<double> pmat{ "problem.mat" ,3,5 };
	for (auto i = 0; i < 3; ++i) {
		pmat.write_row(mat.row(i), i);
	}

	SimplexMethod<double> simp{ "problem.mat",vec_b,vec_c };

	double max_val;
	SimplexMethod<double>::SolutionType sol = move(simp.solve(max_val));
	for (auto iter = sol.begin(); iter != sol.end(); ++iter) {
		cout << "x" << iter->first << " = " << iter->second << "\n";
	}
	cout << "maximum value: " << max_val << "\n";
}

void test_LargeScaleSimplexMethod() {
	// generate a 2000x3000 random matrix
	
	Eigen::MatrixXd* matrix= new Eigen::MatrixXd{ Eigen::MatrixXd::Random(2000,3000)};

	// generate random vec_b and vec_c
	Eigen::MatrixXd vec_b{ Eigen::MatrixXd::Random(2000,1) };
	vec_b = vec_b.cwiseAbs().eval();
	Eigen::MatrixXd vec_c{ Eigen::MatrixXd::Random(1,3000) };

	// write matrix into file
	OnDiskMatrix<double> ondisk{ "matrix.mat",2000,3000 };
	for (auto i = 0; i < matrix->rows(); ++i) {
		ondisk.write_row(matrix->row(i), i);
	}

	// release memory
	delete matrix;

	// run simplex method
	SimplexMethod<double> simplex{ "matrix.mat",vec_b,vec_c };

	double max_val;
	auto sol{ move(simplex.solve(max_val)) };

	cout << max_val;
}

#endif // COMPILE_TEST
