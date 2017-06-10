#ifndef DEF_ONDISKMATRIX_HPP
#define DEF_ONDISKMATRIX_HPP

#include "__include.hpp"
#include "OnDiskMatrixTypeInfo.hpp"


constexpr int type_hint_length = 3;

struct OnDiskMatrixHeader {
	int32_t rows;
	int32_t cols;
	int32_t type_size;
	char type_hint[type_hint_length];
};

constexpr int OnDiskMatrixHeader_size = sizeof(OnDiskMatrixHeader);


template<typename V>
class OnDiskMatrixBase :protected TypeInfo<V> {
public:
	using value_type = V;
	using MatrixType = Eigen::Matrix<value_type, -1, -1,Eigen::RowMajor>;

	// opens a matrix from file
	OnDiskMatrixBase(const string &filename) {
		file.open(filename, ios::binary | ios::in | ios::out);

		// read the header
		file.seekg(ios::beg);
		file.read(reinterpret_cast<char*>(&header), OnDiskMatrixHeader_size);

		// check validity
		verify_type();
	}

	// create a new empty matrix
	OnDiskMatrixBase(const string &filename, int rows, int cols) {
		// initialize header
		header.rows = rows;
		header.cols = cols;
		header.type_size = type_size;
		for (auto i = 0; i < type_hint_length; ++i) {
			header.type_hint[i] = type_hint[i];
		}

		// open file
		file.open(filename, ios::binary|ios::in|ios::out|ios::trunc);
		//file.seekp(ios::beg);

		// write header
		file.write(reinterpret_cast<const char*>(&header), OnDiskMatrixHeader_size);

		// fill the matrix with initialization value
		fill(value_type{});
	}

	virtual ~OnDiskMatrixBase() {
		//file.close();
	};


	// force writing data into file
	void flush() {
		file << std::flush;
	}

	MatrixType read_row(int row_ptr) {
		MatrixType ret{ 1,header.cols };
		file.seekg(get_element_location(row_ptr, 0));
		file.read(reinterpret_cast<char*>(ret.data()), header.cols * type_size);
		return ret;
	}

	// fill the entire matrix with a value
	void fill(const value_type &val) {
		MatrixType row{1,header.cols};
		row.fill(val);
		for (auto i = 0; i < header.rows; ++i) {
			write_row(row, i);
		}
	}

	void write_row(const MatrixType& matrix, int row_ptr) {
		assert(matrix.cols() == header.cols);

		file.seekp(get_element_location(row_ptr, 0));
		file.write(reinterpret_cast<const char*>(matrix.data()), header.cols * type_size);
		flush();
	}

	// avoid using
	value_type get_element(int row_ptr, int col_ptr) {
		value_type ret;
		file.seekg(get_element_location(row_ptr, col_ptr));
		file.read(reinterpret_cast<char*>(&ret), type_size);
		return ret;
	}

	// avoid using
	value_type set_element(const value_type &val, int row_ptr, int col_ptr) {
		file.seekp(get_element_location(row_ptr, col_ptr));
		file.write(reinterpret_cast<const char*>(&val), type_size);
	}

	// generate another matrix, which is the transpose of this matrix
	void generate_transpose_matrix(const string& filename) {
		// generate a new matrix
		OnDiskMatrixBase<value_type> new_matrix{ filename,header.cols,header.rows };
		for (auto i = 0; i < header.rows; ++i) {
			auto row = move(read_row(i));
			row.transposeInPlace();
			new_matrix.write_col(row, i);
		}
	}

	const OnDiskMatrixHeader &get_header() const { return header; }

	int rows() const { return header.rows; }
	int cols() const { return header.cols; }

protected:
	OnDiskMatrixHeader header;
	fstream file;

	streampos get_element_location(int row_ptr, int col_ptr) {
		assert(row_ptr > -1 && row_ptr < header.rows);
		assert(col_ptr > -1 && col_ptr < header.cols);

		return (streampos)OnDiskMatrixHeader_size +
			((streampos)row_ptr * (streampos)header.cols +
			(streampos)col_ptr) * type_size;
	}

	void write_col(const MatrixType& matrix, int col_ptr) {
		assert(matrix.rows() == header.rows);

		for (auto i = 0; i < header.rows; ++i) {
			file.seekp(get_element_location(i, col_ptr));
			file.write(reinterpret_cast<const char*>(&(matrix(i,0))), type_size);
		}
		flush();
	}

	virtual void verify_type() {

		expr_check(header.type_size == type_size, "The size of value type does not match.");
		for (auto i = 0; i < type_hint_length; ++i) {
			expr_check(header.type_hint[i] == type_hint[i], "The type description information does not match.");
		}

	}
};


template<typename V>
class OnDiskMatrix :public OnDiskMatrixBase<V> {
public:
	using base_type = OnDiskMatrixBase<V>;

	OnDiskMatrix(const string &filename) :base_type{ filename } {};
	OnDiskMatrix(const string &filename, int rows, int cols) :base_type{ filename,rows,cols } {}
	OnDiskMatrix(const OnDiskMatrix& other) = delete;
	OnDiskMatrix(OnDiskMatrix&& other) = delete;
	virtual ~OnDiskMatrix(){}
};



#endif // !DEF_ONDISKMATRIX_HPP
