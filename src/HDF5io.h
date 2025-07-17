// HDF5io.h
//
// Breannan Smith, Yonghao Yue
// Last updated: 03/22/2022

#ifndef __HDF5IO_H__
#define __HDF5IO_H__

#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <hdf5.h>
#include <cassert>
#include <vector>

#include <iostream>

inline void tokenize( const std::string& str, const char chr, std::vector<std::string>& tokens )
{
  std::string::size_type substring_start = 0;
  std::string::size_type substring_end = str.find_first_of( chr, substring_start );
  while( substring_end != std::string::npos )
  {
    tokens.emplace_back( str.substr( substring_start, substring_end - substring_start ) );
    substring_start = substring_end + 1;
    substring_end = str.find_first_of( chr, substring_start );
  }
  // Grab the trailing substring, if present
  if( substring_start < str.size() )
  {
    tokens.emplace_back( str.substr( substring_start ) );
  }
  // Case of final character the delimiter
  if( str.back() == chr )
  {
    tokens.emplace_back( "" );
  }
}

inline std::vector<std::string> tokenize( const std::string& str, const char delimiter )
{
  std::vector<std::string> tokens;
  tokenize( str, delimiter, tokens );
  return tokens;
}

enum class HDF5AccessType : std::uint8_t
{
  READ_ONLY,
  READ_WRITE
};

namespace HDF5SupportedTypes
{
  //see http://davis.lbl.gov/Manuals/HDF5-1.8.7/UG/11_Datatypes.html
  template <typename T> constexpr bool isSupportedEigenType() { return false; }
  template <> constexpr bool isSupportedEigenType<char>() { return true; }
  template <> constexpr bool isSupportedEigenType<signed char>() { return true; }
  template <> constexpr bool isSupportedEigenType<unsigned char>() { return true; }
  template <> constexpr bool isSupportedEigenType<short>() { return true; }
  template <> constexpr bool isSupportedEigenType<unsigned short>() { return true; }
  template <> constexpr bool isSupportedEigenType<int>() { return true; }
  template <> constexpr bool isSupportedEigenType<unsigned>() { return true; }
  template <> constexpr bool isSupportedEigenType<long>() { return true; }
  template <> constexpr bool isSupportedEigenType<unsigned long>() { return true; }
  template <> constexpr bool isSupportedEigenType<long long>() { return true; }
  template <> constexpr bool isSupportedEigenType<unsigned long long>() { return true; }
  template <> constexpr bool isSupportedEigenType<float>() { return true; }
  template <> constexpr bool isSupportedEigenType<double>() { return true; }
  template <> constexpr bool isSupportedEigenType<long double>() { return true; }
  //template <> constexpr bool isSupportedEigenType<hsize_t>() { return true; } == unsigned long long
  //template <> constexpr bool isSupportedEigenType<hssize_t>() { return true; } == long long
  //template <> constexpr bool isSupportedEigenType<herr_t>() { return true; } == int
  template <> constexpr bool isSupportedEigenType<hbool_t>() { return true; }
}

// Wrapper for HDF5 group id. Functionality for HDF5 types differs only in the required close operaiton.
template <herr_t H5CloseOperation( hid_t id )>
class HDFID final
{
public:
  HDFID(): m_hid_t( -1 ) {}
  explicit HDFID( const hid_t value ): m_hid_t( value ) {}
  HDFID( HDFID&& other ): m_hid_t( other.m_hid_t ) { other.m_hid_t = -1; }
	~HDFID(){ if( m_hid_t >= 0 ) H5CloseOperation( m_hid_t ); }

  HDFID& operator=( HDFID&& other )
  {
    // NB: Backup hid_t to guard against self assignment
    const hid_t others_hid_t{ other.m_hid_t };
    other.m_hid_t = -1;
    m_hid_t = others_hid_t;
    return *this;
  }

  operator hid_t() const
  {
    return m_hid_t;
  }

private:
  HDFID( const HDFID& ) = delete;
  HDFID& operator=( const HDFID& ) = delete;

  hid_t m_hid_t;
};

using HDFGID = HDFID<H5Gclose>;
using HDFSID = HDFID<H5Sclose>;
using HDFDID = HDFID<H5Dclose>;
using HDFTID = HDFID<H5Tclose>;

class HDF5File final
{
	HDF5File() = delete;
public:
	HDF5File( const std::string& file_name, const HDF5AccessType& access_type )
  : m_hdf_file_id( -1 )
  , m_file_opened( false )
	{
		// Attempt to open a file
	  switch( access_type )
	  {
	    case HDF5AccessType::READ_WRITE:
	      m_hdf_file_id = H5Fcreate( file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
	      break;
	    case HDF5AccessType::READ_ONLY:
	      m_hdf_file_id = H5Fopen( file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	      break;
	  }
	  // Check that the file successfully opened
	  if( m_hdf_file_id < 0 )
	  {
	    throw std::string{ "Failed to open file: " } + file_name;
	  }
	  m_file_opened = true;
	}

  ~HDF5File()
	{
		if( m_file_opened )
	  {
	    assert( m_hdf_file_id >= 0 );
	    H5Fclose( m_hdf_file_id );
	  }
	}

	bool isOpened()
	{
		return m_file_opened;
	}
  
  bool doesGroupExist( const std::string& group_name ) const
  {
    return H5Lexists( m_hdf_file_id, group_name.c_str(), H5P_DEFAULT );
  }

  bool doesDatasetExist( const std::string& group_name, const std::string& dataset_name ) const
  {
    const HDFGID grp_id{ findGroup( group_name ) };
    return H5Lexists( grp_id, dataset_name.c_str(), H5P_DEFAULT );
  }

  template <typename Scalar>
	void writeScalar( const std::string& group_name, const std::string& variable_name, const Scalar& variable ) const
  {
    //Eigen::Matrix<double,1,1> output_mat;
    //output_mat << variable;
    //writeMatrix( group_name, variable_name, output_mat );
    const HDFSID dataspace_id{ H5Screate( H5S_SCALAR ) };
    if( dataspace_id < 0 )
    {
      throw std::string{ "Failed to create HDF data space" };
    }

    // Open the requested group
    const HDFGID grp_id{ findOrCreateGroup( group_name ) };

    const HDFDID dataset_id{ H5Dcreate2( grp_id, variable_name.c_str(), computeHDFType<Scalar>(), dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ) };
    if( dataset_id < 0 )
    {
      throw std::string{ "Failed to create HDF data set" };
    }

    const herr_t status_write{ H5Dwrite( dataset_id, computeHDFType<Scalar>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &variable ) };
    if( status_write < 0 )
    {
      throw std::string{ "Failed to write HDF data" };
    }
  }

  template <typename Scalar>
  void readScalar( const std::string& group_name, const std::string& variable_name, Scalar& variable ) const
  {
    // Open the requested group
    const HDFGID grp_id{ findGroup( group_name ) };

    const HDFDID dataset_id{ H5Dopen2( grp_id, variable_name.c_str(), H5P_DEFAULT ) };
    if( dataset_id < 0 )
    {
      throw std::string{ "Failed to open HDF data set" };
    }
    //if( getNativeType( dataset_id ) != computeHDFType<Scalar>() )
    if( !H5Tequal( getNativeType( dataset_id ), computeHDFType<Scalar>() ) )
    {
      throw std::string{ "Requested HDF data set is not of given type from Eigen variable" };
    }

    Eigen::ArrayXi dimensions;
    getDimensions( dataset_id, dimensions );
    if( dimensions.size() == 2 )
    {
      Eigen::Matrix<Scalar, 1, 1> input_mat;
      readMatrix( group_name, variable_name, input_mat );
      variable = input_mat( 0, 0 );
      return;
    }
    else if( dimensions.size() != 0 )
    {
      throw std::string{ "Invalid dimensions for scalar type in file" };
    }

    const herr_t read_status{ H5Dread( dataset_id, getNativeType( dataset_id ), H5S_ALL, H5S_ALL, H5P_DEFAULT, &variable ) };
    if( read_status < 0 )
    {
      throw std::string{ "Failed to read data from HDF file" };
    }
  }

	template <typename Derived>
  void writeMatrix( const std::string& group_name, const std::string& variable_name, const Eigen::DenseBase<Derived>& eigen_variable ) const
  {
    static_assert( HDF5SupportedTypes::isSupportedEigenType<typename Derived::Scalar>(), "Error, scalar type of Eigen variable must be float, double, uint8_t, uint32_t, uint64_t, int8_t, int32_t, or int64_t" );

    assert( eigen_variable.rows() >= 0 ); assert( eigen_variable.cols() >= 0 );
    const hsize_t dims[2] = { hsize_t( eigen_variable.rows() ), hsize_t( eigen_variable.cols() ) };
    const HDFSID dataspace_id{ H5Screate_simple( 2, dims, nullptr ) };
    if( dataspace_id < 0 )
    {
      throw std::string{ "Failed to create HDF data space" };
    }

    // Open the requested group
    const HDFGID grp_id{ findOrCreateGroup( group_name ) };

    const HDFDID dataset_id{ H5Dcreate2( grp_id, variable_name.c_str(), computeHDFType( eigen_variable ), dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ) };
    if( dataset_id < 0 )
    {
      throw std::string{ "Failed to create HDF data set" };
    }

    // Convert to row major format, if needed
    Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> col_major_output_data;
    if( isColumnMajor( eigen_variable ) )
    {
      col_major_output_data.resize( eigen_variable.rows(), eigen_variable.cols() );
      col_major_output_data = eigen_variable.derived().matrix();
    }

    const herr_t status_write{ H5Dwrite( dataset_id, computeHDFType( eigen_variable ), H5S_ALL, H5S_ALL, H5P_DEFAULT, isColumnMajor( eigen_variable ) ? col_major_output_data.data() : eigen_variable.derived().data() ) };
    if( status_write < 0 )
    {
      throw std::string{ "Failed to write HDF data" };
    }
  }

  
  /*
   * python対応用の出力形式
   * 2*2行列の場合： 第１列と第２列を行列の総数分だけ開けて追加
   * (例: index: 0, 行列総数 10, 行列 A を追加  ->  追加場所 eigen_variable.col(0) = A.col(0), eigen_variable.col(10) = A.col(1))
   * 第１列と第２列を並べて追加した場合正しい出力にならないため要注意
  */
  template <typename Derived>
  void writeSquareMatrixArray( const std::string& group_name, const std::string& variable_name, const Eigen::DenseBase<Derived>& eigen_variable) const
  {
    static_assert( HDF5SupportedTypes::isSupportedEigenType<typename Derived::Scalar>(), "Error, scalar type of Eigen variable must be float, double, uint8_t, uint32_t, uint64_t, int8_t, int32_t, or int64_t" );

    assert( eigen_variable.rows() >= 0 ); assert( eigen_variable.cols() >= 0 ); assert( eigen_variable.cols() % eigen_variable.rows() == 0.0 );
    const hsize_t dims[3] = { hsize_t( eigen_variable.rows() ), hsize_t( eigen_variable.rows() ), hsize_t( eigen_variable.cols() / eigen_variable.rows() ) };
    const HDFSID dataspace_id{ H5Screate_simple( 3, dims, nullptr ) };
    if( dataspace_id < 0 )
    {
      throw std::string{ "Failed to create HDF data space" };
    }

    // Open the requested group
    const HDFGID grp_id{ findOrCreateGroup( group_name ) };

    const HDFDID dataset_id{ H5Dcreate2( grp_id, variable_name.c_str(), computeHDFType( eigen_variable ), dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ) };
    if( dataset_id < 0 )
    {
      throw std::string{ "Failed to create HDF data set" };
    }

    // Convert to row major format, if needed
    Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> col_major_output_data;
    if( isColumnMajor( eigen_variable ) )
    {
      col_major_output_data.resize( eigen_variable.rows(), eigen_variable.cols() );
      col_major_output_data = eigen_variable.derived().matrix();
    }

    const herr_t status_write{ H5Dwrite( dataset_id, computeHDFType( eigen_variable ), H5S_ALL, H5S_ALL, H5P_DEFAULT, isColumnMajor( eigen_variable ) ? col_major_output_data.data() : eigen_variable.derived().data() ) };
    if( status_write < 0 )
    {
      throw std::string{ "Failed to write HDF data" };
    }
  }
  
  template <typename Derived>
  void writeVector( const std::string& group_name, const std::string& variable_name, const Eigen::DenseBase<Derived>& eigen_variable ) const
  {
    static_assert( HDF5SupportedTypes::isSupportedEigenType<typename Derived::Scalar>(), "Error, scalar type of Eigen variable must be float, double, uint8_t, uint32_t, uint64_t, int8_t, int32_t, or int64_t" );

    assert( eigen_variable.rows() >= 0 ); assert( eigen_variable.cols() == 1 );
    const hsize_t dims[1] = { hsize_t( eigen_variable.rows() ) };
    const HDFSID dataspace_id{ H5Screate_simple( 1, dims, nullptr ) };
    if( dataspace_id < 0 )
    {
      throw std::string{ "Failed to create HDF data space" };
    }

    // Open the requested group
    const HDFGID grp_id{ findOrCreateGroup( group_name ) };

    const HDFDID dataset_id{ H5Dcreate2( grp_id, variable_name.c_str(), computeHDFType( eigen_variable ), dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ) };
    if( dataset_id < 0 )
    {
      throw std::string{ "Failed to create HDF data set" };
    }

    // Convert to row major format, if needed

    const herr_t status_write{ H5Dwrite( dataset_id, computeHDFType( eigen_variable ), H5S_ALL, H5S_ALL, H5P_DEFAULT, eigen_variable.derived().data() ) };
    if( status_write < 0 )
    {
      throw std::string{ "Failed to write HDF data" };
    }
  }

  template <typename Derived>
  void readMatrix( const std::string& group_name, const std::string& variable_name, Eigen::DenseBase<Derived>& eigen_variable ) const
  {
    static_assert( HDF5SupportedTypes::isSupportedEigenType<typename Derived::Scalar>(), "Error, scalar type of Eigen variable must be float, double, uint8_t, uint32_t, uint64_t, int8_t, int32_t, or int64_t" );

    // Open the requested group
    const HDFGID grp_id{ findGroup( group_name ) };

    const HDFDID dataset_id{ H5Dopen2( grp_id, variable_name.c_str(), H5P_DEFAULT ) };
    if( dataset_id < 0 )
    {
      throw std::string{ "Failed to open HDF data set" };
    }

    //std::cout << "types: " << nativeTypeString( getNativeType( dataset_id ) ) << ", " << nativeTypeString( computeHDFType( eigen_variable ) ) << std::endl;

    if( !H5Tequal( getNativeType( dataset_id ), computeHDFType( eigen_variable ) ) )
    {
      throw std::string{ "Requested HDF data set is not of given type from Eigen variable" };
    }

    Eigen::ArrayXi dimensions;
    getDimensions( dataset_id, dimensions );
    if( dimensions.size() == 0 )
    {
      throw std::string{ "The data is a scalar, not a matrix or vector. Try readdouble() instead." };
    }

    if( ( dimensions < 0 ).any() )
    {
      throw std::string{ "Negative dimensions for Eigen matrix type in file" };
    }

    if( dimensions.size() == 1 )
    {
      if( eigen_variable.cols() != 1 )
      {
        throw std::string{ "Eigen type is not a vector" };
      }

      if( rowsFixed( eigen_variable ) && eigen_variable.rows() != dimensions( 0 ) )
      {
        throw std::string{ "Eigen type of fixed row size does not have correct number of rows" };
      }

      // Resize the Eigen type
      eigen_variable.derived().resize( dimensions( 0 ), 1 );

      const herr_t read_status{ H5Dread( dataset_id, getNativeType( dataset_id ), H5S_ALL, H5S_ALL, H5P_DEFAULT, eigen_variable.derived().data() ) };
      if( read_status < 0 )
      {
        throw std::string{ "Failed to read data from HDF file" };
      }
    }
    else if( dimensions.size() == 2 )
    {
      // If the Eigen type has a fixed dimension, ensure it is correct
      if( rowsFixed( eigen_variable ) && eigen_variable.rows() != dimensions( 0 ) )
      {
        throw std::string{ "Eigen type of fixed row size does not have correct number of rows" };
      }
      if( colsFixed( eigen_variable ) && eigen_variable.cols() != dimensions( 1 ) )
      {
        throw std::string{ "Eigen type of fixed cols size does not have correct number of cols" };
      }
      // Resize the Eigen type
      eigen_variable.derived().resize( dimensions( 0 ), dimensions( 1 ) );

      // Convert to row major format, if needed
      Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> col_major_input_data;
      if( isColumnMajor( eigen_variable ) )
      {
        col_major_input_data.resize( eigen_variable.rows(), eigen_variable.cols() );
      }

      const herr_t read_status{ H5Dread( dataset_id, getNativeType( dataset_id ), H5S_ALL, H5S_ALL, H5P_DEFAULT, isColumnMajor( eigen_variable ) ? col_major_input_data.data() : eigen_variable.derived().data() ) };
      if( read_status < 0 )
      {
        throw std::string{ "Failed to read data from HDF file" };
      }
      if( isColumnMajor( eigen_variable ) )
      {
        eigen_variable.derived().matrix() = col_major_input_data;
      }
    }
    else
    {
      throw std::string{ "Invalid dimensions for Eigen matrix type in file" };
    }
  }

  template <typename Derived>
  void writeSparseMatrix( const std::string& group_name, const std::string& variable_name, const Eigen::SparseMatrixBase<Derived>& sparse_matrix ) const
  {
    static_assert( !Derived::IsRowMajor, "Error, HDF5 sparse matrix output only supported for column major matrices" );

    // Grab the sparse data
    Eigen::Matrix<typename Derived::StorageIndex,Eigen::Dynamic,1> col_ptr;
    Eigen::Matrix<typename Derived::StorageIndex,Eigen::Dynamic,1> row_ind;
    Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,1> val;
    extractDataCCS( sparse_matrix, col_ptr, row_ind, val );

    // Write out the sparse data
    const Eigen::Index nrows{ sparse_matrix.rows() };
    writeScalar( group_name, variable_name + "_nrows", nrows );
    const Eigen::Index ncols{ sparse_matrix.cols() };
    writeScalar( group_name, variable_name + "_ncols", ncols );
    writeMatrix( group_name, variable_name + "_col_ptr", col_ptr );
    writeMatrix( group_name, variable_name + "_row_ind", row_ind );
    writeMatrix( group_name, variable_name + "_val", val );
  }

  template <typename Derived>
  void readSparseMatrix( const std::string& group_name, const std::string& variable_name, Eigen::SparseMatrixBase<Derived>& sparse_matrix ) const
  {
    static_assert( !Derived::IsRowMajor, "Error, HDF5 sparse matrix input only supported for column major matrices" );

    // Read in the sparse data
    std::vector< Eigen::Triplet<typename Derived::Scalar> > triplet_list;
    {
      typename Derived::Index nrows;
      readScalar( group_name, variable_name + "_nrows", nrows );
      typename Derived::Index ncols;
      readScalar( group_name, variable_name + "_ncols", ncols );
      Eigen::Matrix<typename Derived::Index,Eigen::Dynamic,1> col_ptr;
      readMatrix( group_name, variable_name + "_col_ptr", col_ptr );
      if( col_ptr.size() != ncols + 1 )
      {
        throw std::string{ "Column count and size of col_ptr do not agree for sparse matrix" };
      }
      Eigen::Matrix<typename Derived::Index,Eigen::Dynamic,1> row_ind;
      readMatrix( group_name, variable_name + "_row_ind", row_ind );
      Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,1> val;
      readMatrix( group_name, variable_name + "_val", val );
      if( row_ind.size() != val.size() )
      {
        throw std::string{ "Row indices and value size do not agree for sparse matrix" };
      }
      sparse_matrix.derived().resize( nrows, ncols );
      sparse_matrix.derived().reserve( int( val.size() ) );
      for( typename Derived::Index col_num = 0; col_num < ncols; ++col_num )
      {
        assert( col_ptr( col_num ) <= col_ptr( col_num + 1 ) );
        for( typename Derived::Index idx = col_ptr( col_num ); idx < col_ptr( col_num + 1 ); ++idx )
        {
          typename Derived::Index row_num = row_ind( idx );
          triplet_list.push_back( Eigen::Triplet<typename Derived::Scalar>( row_num, col_num, val( idx ) ) );
        }
      }
    }
    sparse_matrix.derived().setFromTriplets( triplet_list.cbegin(), triplet_list.cend() );
    sparse_matrix.derived().makeCompressed();
  }

  template <typename Derived>
  void readSquareMatrixArray( const std::string& group_name, const std::string& variable_name, Eigen::DenseBase<Derived>& eigen_variable ) const
  {
    static_assert( HDF5SupportedTypes::isSupportedEigenType<typename Derived::Scalar>(), "Error, scalar type of Eigen variable must be float, double, uint8_t, uint32_t, uint64_t, int8_t, int32_t, or int64_t" );

    // Open the requested group
    const HDFGID grp_id{ findGroup( group_name ) };

    const HDFDID dataset_id{ H5Dopen2( grp_id, variable_name.c_str(), H5P_DEFAULT ) };
    if( dataset_id < 0 )
    {
      throw std::string{ "Failed to open HDF data set" };
    }

    //std::cout << "types: " << nativeTypeString( getNativeType( dataset_id ) ) << ", " << nativeTypeString( computeHDFType( eigen_variable ) ) << std::endl;

    if( !H5Tequal( getNativeType( dataset_id ), computeHDFType( eigen_variable ) ) )
    {
      throw std::string{ "Requested HDF data set is not of given type from Eigen variable" };
    }

    Eigen::ArrayXi dimensions;
    getDimensions( dataset_id, dimensions );
    
    if( dimensions.size() == 0 )
    {
      throw std::string{ "The data is a scalar, not a matrix or vector. Try readdouble() instead." };
    }

    if( ( dimensions < 0 ).any() )
    {
      throw std::string{ "Negative dimensions for Eigen matrix type in file" };
    }

    if( dimensions.size() == 3 )
    {
      // If the Eigen type has a fixed dimension, ensure it is correct
      if( rowsFixed( eigen_variable ) && eigen_variable.rows() != dimensions( 0 ) )
      {
        throw std::string{ "Eigen type of fixed row size does not have correct number of rows" };
      }
      if( colsFixed( eigen_variable ) && eigen_variable.cols() != dimensions( 1 ) )
      {
        throw std::string{ "Eigen type of fixed cols size does not have correct number of cols" };
      }
      // Resize the Eigen type
      eigen_variable.derived().resize( dimensions( 0 ), dimensions( 1 ) * dimensions( 2 ) );

      // Convert to row major format, if needed
      Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> col_major_input_data;
      if( isColumnMajor( eigen_variable ) )
      {
        col_major_input_data.resize( eigen_variable.rows(), eigen_variable.cols() );
      }

      const herr_t read_status{ H5Dread( dataset_id, getNativeType( dataset_id ), H5S_ALL, H5S_ALL, H5P_DEFAULT, isColumnMajor( eigen_variable ) ? col_major_input_data.data() : eigen_variable.derived().data() ) };
      if( read_status < 0 )
      {
        throw std::string{ "Failed to read data from HDF file" };
      }
      if( isColumnMajor( eigen_variable ) )
      {
        eigen_variable.derived().matrix() = col_major_input_data;
      }
    }
    else
    {
      throw std::string{ "Invalid dimensions for Eigen matrix type in file" };
    }
  }
  
  void writeString( const std::string& group_name, const std::string& variable_name, const std::string& string_variable ) const
  {
    //const char* string_data[1] = { string_variable.c_str() };
    int str_len = string_variable.length();
    //const HDFTID file_type{ H5Tcopy( H5T_STRING ) };
    const HDFTID file_type{ H5Tcopy( H5T_C_S1 ) };
    if( file_type < 0 )
    {
      throw std::string{ "Failed to create HDF5 string file type" };
    }

    const herr_t set_file_pad_status{ H5Tset_strpad( file_type, H5T_STR_NULLPAD ) };
    if( set_file_pad_status < 0 )
    {
      throw std::string{ "Failed to set HDF5 string file pad" };
    }

    const herr_t set_file_cset_status{ H5Tset_cset( file_type, H5T_CSET_ASCII ) };
    if( set_file_cset_status < 0 )
    {
      throw std::string{ "Failed to set HDF5 string file cset" };
    }

    const herr_t set_file_size_status{ H5Tset_size( file_type, str_len ) };
    if( set_file_size_status < 0 )
    {
      throw std::string{ "Failed to set HDF5 string file size" };
    }

    const hsize_t dims[1] = { 1 };
    const HDFSID space{ H5Screate_simple( 1, dims, nullptr ) };
    if( space < 0 )
    {
      throw std::string{ "Failed to create HDF space" };
    }

    // Open the requested group
    const HDFGID grp_id{ findOrCreateGroup( group_name ) };

    const HDFDID dset{ H5Dcreate2( grp_id, variable_name.c_str(), file_type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ) };
    if( dset < 0 )
    {
      throw std::string{ "Failed to create HDF dataset" };
    }
    const herr_t write_status{ H5Dwrite( dset, file_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, string_variable.c_str() ) };
    if( write_status < 0 )
    {
      throw std::string{ "Failed to write HDF data" };
    }
  }

  void readString( const std::string& group_name, const std::string& variable_name, std::string& string_variable ) const
  {
    const HDFGID grp_id{ findGroup( group_name ) };

    const HDFDID dset{ H5Dopen2( grp_id, variable_name.c_str(), H5P_DEFAULT ) };
    if( dset < 0 )
    {
      throw std::string{ "Failed to open HDF dataset" };
    }
    const HDFTID file_type{ H5Dget_type( dset ) };
    if( file_type < 0 )
    {
      throw std::string{ "Failed to get HDF file type" };
    }

    H5T_class_t type_class = H5Tget_class( file_type );
    if( type_class != H5T_STRING )
    {
      throw std::string{ "Dataset class is not H5T_STRING" };
    }

    hid_t cset = H5Tget_cset( file_type );
    if( cset != H5T_CSET_ASCII )
    {
      throw std::string{ "Currently only support ASCII string" };
    }

    const HDFTID native_type{ ( type_class == H5T_BITFIELD ) ? H5Tcopy( file_type ) : H5Tget_native_type( file_type, H5T_DIR_DEFAULT ) };
    size_t n_size = std::max( H5Tget_size( file_type ), H5Tget_size( native_type ) );

    char* data = (char*)calloc( n_size + 1, sizeof(char) );
    if( data == nullptr )
    {
      throw std::string{ "Unable to allocate memory" };
    }

    const HDFTID mem_type{ H5Tcopy( H5T_C_S1 ) };
    H5Tset_cset( mem_type, H5T_CSET_ASCII );
    H5Tset_size( mem_type, n_size + 1 );

    const herr_t read_status{ H5Dread( dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data ) };
    if( read_status < 0 )
    {
      free( data );
      throw std::string{ "Failed to read HDF data" };
    }

    string_variable = data;
    free( data );
  }

	void enumerateGroups( const std::string& group_name, std::vector<std::string>& out_group_names ) const
	{
		const HDFGID grp_id{ findGroup( group_name ) };
		out_group_names.clear();
		H5Literate( hid_t( grp_id ), H5_INDEX_NAME, H5_ITER_NATIVE, NULL, iterate_group_func, (void *)&out_group_names );
	}

	void enumerateDataSets( const std::string& group_name, std::vector<std::string>& out_dataset_names ) const
	{
		const HDFGID grp_id{ findGroup( group_name ) };
		out_dataset_names.clear();
		H5Literate( hid_t( grp_id ), H5_INDEX_NAME, H5_ITER_NATIVE, NULL, iterate_dataset_func, (void *)&out_dataset_names );
	}

private:
	static herr_t iterate_group_func( hid_t loc_id, const char *name, const H5L_info_t *info,
            void *operator_data )
	{
		//https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5G/h5ex_g_traverse.c
		std::vector<std::string>* out_data = reinterpret_cast<std::vector<std::string>*>( operator_data );

		H5O_info_t infobuf;
		H5Oget_info_by_name( loc_id, name, &infobuf, H5O_INFO_BASIC, H5P_DEFAULT );
		if( infobuf.type == H5O_TYPE_GROUP )
		{
			out_data->push_back( std::string( name ) );
		}

		return 0;
	}

	static herr_t iterate_dataset_func( hid_t loc_id, const char *name, const H5L_info_t *info,
            void *operator_data )
	{
		//https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5G/h5ex_g_traverse.c
		std::vector<std::string>* out_data = reinterpret_cast<std::vector<std::string>*>( operator_data );

		H5O_info_t infobuf;
		H5Oget_info_by_name( loc_id, name, &infobuf, H5O_INFO_BASIC, H5P_DEFAULT );
		if( infobuf.type == H5O_TYPE_DATASET )
		{
			out_data->push_back( std::string( name ) );
		}

		return 0;
	}

	static hid_t getNativeType( const hid_t dataset_id )
  {
    // see http://davis.lbl.gov/Manuals/HDF5-1.8.7/UG/11_Datatypes.html
    const HDFTID dtype_id{ H5Dget_type( dataset_id ) };
    if( dtype_id < 0 )
    {
      throw std::string{ "Failed to get HDF5 datatype from dataset" };
    }
    if( H5Tequal( dtype_id, H5T_NATIVE_CHAR ) > 0 )
    {
      return H5T_NATIVE_CHAR;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_SCHAR ) > 0 )
    {
      return H5T_NATIVE_SCHAR;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_UCHAR ) > 0 )
    {
      return H5T_NATIVE_UCHAR;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_SHORT ) > 0 )
    {
      return H5T_NATIVE_SHORT;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_USHORT ) > 0 )
    {
      return H5T_NATIVE_USHORT;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_INT ) > 0 )
    {
      return H5T_NATIVE_INT;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_UINT ) > 0 )
    {
      return H5T_NATIVE_UINT;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_LONG ) > 0 )
    {
      return H5T_NATIVE_LONG;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_ULONG ) > 0 )
    {
      return H5T_NATIVE_ULONG;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_LLONG ) > 0 )
    {
      return H5T_NATIVE_LLONG;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_ULLONG ) > 0 )
    {
      return H5T_NATIVE_ULLONG;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_FLOAT ) > 0 )
    {
      return H5T_NATIVE_FLOAT;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_DOUBLE ) > 0 )
    {
      return H5T_NATIVE_DOUBLE;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_LDOUBLE ) > 0 )
    {
      return H5T_NATIVE_LDOUBLE;
    }
    /*
    else if( H5Tequal( dtype_id, H5T_NATIVE_HSIZE ) > 0 ) == unsigned long long
    {
      return H5T_NATIVE_HSIZE;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_HSSIZE ) > 0 ) == long long
    {
      return H5T_NATIVE_HSSIZE;
    }
    else if( H5Tequal( dtype_id, H5T_NATIVE_HERR ) > 0 ) == int
    {
      return H5T_NATIVE_HERR;
    }
    //*/
    else if( H5Tequal( dtype_id, H5T_NATIVE_HBOOL ) > 0 )
    {
      return H5T_NATIVE_HBOOL;
    }
    else return -1;
  }

  static std::string nativeTypeString( const hid_t native_type )
  {
    if( H5Tequal( native_type, H5T_NATIVE_CHAR ) > 0 )
    {
      return std::string{ "H5T_NATIVE_CHAR" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_SCHAR ) > 0 )
    {
      return std::string{ "H5T_NATIVE_SCHAR" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_UCHAR ) > 0 )
    {
      return std::string{ "H5T_NATIVE_UCHAR" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_SHORT ) > 0 )
    {
      return std::string{ "H5T_NATIVE_SHORT" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_USHORT ) > 0 )
    {
      return std::string{ "H5T_NATIVE_USHORT" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_INT ) > 0 )
    {
      return std::string{ "H5T_NATIVE_INT" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_UINT ) > 0 )
    {
      return std::string{ "H5T_NATIVE_UINT" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_LONG ) > 0 )
    {
      return std::string{ "H5T_NATIVE_LONG" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_ULONG ) > 0 )
    {
      return std::string{ "H5T_NATIVE_ULONG" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_LLONG ) > 0 )
    {
      return std::string{ "H5T_NATIVE_LLONG" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_ULLONG ) > 0 )
    {
      return std::string{ "H5T_NATIVE_ULLONG" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_FLOAT ) > 0 )
    {
      return std::string{ "H5T_NATIVE_FLOAT" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_DOUBLE ) > 0 )
    {
      return std::string{ "H5T_NATIVE_DOUBLE" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_LDOUBLE ) > 0 )
    {
      return std::string{ "H5T_NATIVE_LDOUBLE" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_HSIZE ) > 0 )
    {
      return std::string{ "H5T_NATIVE_HSIZE" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_HSSIZE ) > 0 )
    {
      return std::string{ "H5T_NATIVE_HSSIZE" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_HERR ) > 0 )
    {
      return std::string{ "H5T_NATIVE_HERR" };
    }
    else if( H5Tequal( native_type, H5T_NATIVE_HBOOL ) > 0 )
    {
      return std::string{ "H5T_NATIVE_HBOOL" };
    }
    else return std::string{ "-1" };
  }

  static void getDimensions( const hid_t dataset_id, Eigen::ArrayXi& dimensions )
  {
    const HDFSID space_id{ H5Dget_space( dataset_id ) };
    if( space_id < 0 )
    {
      throw std::string{ "Failed to open data space" };
    }
    const int rank{ H5Sget_simple_extent_ndims( space_id ) };
    if( rank < 0 )
    {
      throw std::string{ "Failed to get rank" };
    }
    dimensions.resize( rank );
    std::vector<hsize_t> dims( static_cast<std::vector<hsize_t>::size_type>( rank ) );
    assert( int( dims.size() ) == rank );
    const herr_t status_get_simple_extent_dims{ H5Sget_simple_extent_dims( space_id, dims.data(), nullptr ) };
    if( status_get_simple_extent_dims < 0 )
    {
      throw std::string{ "Failed to get extents" };
    }
    for( int i = 0; i < rank; ++i )
    {
      dimensions( i ) = int( dims[i] );
    }
  }

	template <typename ScalarType>
  static constexpr hid_t computeHDFType()
  {
    //see http://davis.lbl.gov/Manuals/HDF5-1.8.7/UG/11_Datatypes.html
    return std::is_same<ScalarType, char>::value ? H5T_NATIVE_CHAR :
      std::is_same<ScalarType, signed char>::value ? H5T_NATIVE_SCHAR :
      std::is_same<ScalarType, unsigned char>::value ? H5T_NATIVE_UCHAR :
      std::is_same<ScalarType, short>::value ? H5T_NATIVE_SHORT :
      std::is_same<ScalarType, unsigned short>::value ? H5T_NATIVE_USHORT :
      std::is_same<ScalarType, int>::value ? H5T_NATIVE_INT :
      std::is_same<ScalarType, unsigned>::value ? H5T_NATIVE_UINT :
      std::is_same<ScalarType, long>::value ? H5T_NATIVE_LONG :
      std::is_same<ScalarType, unsigned long>::value ? H5T_NATIVE_ULONG :
      std::is_same<ScalarType, long long>::value ? H5T_NATIVE_LLONG :
      std::is_same<ScalarType, unsigned long long>::value ? H5T_NATIVE_ULLONG :
      std::is_same<ScalarType, float>::value ? H5T_NATIVE_FLOAT :
      std::is_same<ScalarType, double>::value ? H5T_NATIVE_DOUBLE :
      std::is_same<ScalarType, long double>::value ? H5T_NATIVE_LDOUBLE :
      //std::is_same<ScalarType, hsize_t>::value ? H5T_NATIVE_HSIZE : == unsigned long long
      //std::is_same<ScalarType, hssize_t>::value ? H5T_NATIVE_HSSIZE : == long long
      //std::is_same<ScalarType, herr_t>::value ? H5T_NATIVE_HERR : == int
      std::is_same<ScalarType, hbool_t>::value ? H5T_NATIVE_HBOOL :
      -1;
  }

	template <typename Derived>
  static constexpr hid_t computeHDFType( const Eigen::EigenBase<Derived>& )
  {
    return computeHDFType<typename Derived::Scalar>();
  }

  template <typename Derived>
  static constexpr bool isColumnMajor( const Eigen::EigenBase<Derived>& )
  {
    return !Derived::IsRowMajor;
  }

  template <typename Derived>
  static constexpr bool rowsFixed( const Eigen::EigenBase<Derived>& )
  {
    return Derived::RowsAtCompileTime != Eigen::Dynamic;
  }

  template <typename Derived>
  static constexpr bool colsFixed( const Eigen::EigenBase<Derived>& )
  {
    return Derived::ColsAtCompileTime != Eigen::Dynamic;
  }

	HDFID<H5Gclose> findGroup( const std::string& group_name ) const
	{
		HDFGID group_id;
	  if( group_name.empty() )
	  {
	    group_id = HDFGID{ H5Gopen2( m_hdf_file_id, "/", H5P_DEFAULT ) };
	  }
	  else
	  {
	    group_id = HDFGID{ H5Gopen2( m_hdf_file_id, group_name.c_str(), H5P_DEFAULT ) };
	  }
	  if( group_id < 0 )
	  {
	    throw std::string{ "Failed to find group: " } + group_name;
	  }
	  return group_id;
	}

	HDFGID findOrCreateGroup( const std::string& group_name ) const
	{
	  if( group_name.empty() )
	  {
	    return HDFGID{ H5Gopen2( m_hdf_file_id, "/", H5P_DEFAULT ) };
	  }

	  const std::vector<std::string> components{ tokenize( group_name, '/' ) };

	  // For each level in the path
	  std::string current_depth;
	  for (const std::string& val : components )
	  {
	    current_depth = current_depth + '/' + val;
	    // If the level does not exist
	    if( 0 == H5Lexists( m_hdf_file_id, current_depth.c_str(), H5P_DEFAULT ) )
	    {
	      // Create the level
	      HDFGID group_id = HDFGID{ H5Gcreate2( m_hdf_file_id, current_depth.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ) };
	      if( group_id < 0 )
	      {
	        throw std::string{ "Failed to create group: " } + current_depth;
	      }
	    }
	  }

	  // Return the final constructed level
	  HDFGID group_id = HDFGID{ H5Gopen2( m_hdf_file_id, group_name.c_str(), H5P_DEFAULT ) };
	  if( group_id < 0 )
	  {
	    throw std::string{ "Failed to open group: " } + group_name;
	  }

	  return group_id;
	}

  // TODO: Don't extract data first, just directly save out to disk
  template<typename Derived>
  static void extractDataCCS( const Eigen::SparseMatrixBase<Derived>& A, Eigen::Matrix<typename Derived::StorageIndex,Eigen::Dynamic,1>& outer_ptr, Eigen::Matrix<typename Derived::StorageIndex,Eigen::Dynamic,1>& inner_ptr, Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,1>& val )
  {
    outer_ptr = Eigen::Map< const Eigen::Matrix<typename Derived::StorageIndex,Eigen::Dynamic,1> >( A.derived().outerIndexPtr(), A.outerSize() + 1 );
    inner_ptr = Eigen::Map< const Eigen::Matrix<typename Derived::StorageIndex,Eigen::Dynamic,1> >( A.derived().innerIndexPtr(), A.derived().nonZeros() );
    val = Eigen::Map< const Eigen::Array<typename Derived::Scalar,Eigen::Dynamic,1> >( A.derived().valuePtr(), A.derived().nonZeros() );
  }

  hid_t fileID() const
  {
    return m_hdf_file_id;
  }

	hid_t m_hdf_file_id;
  bool m_file_opened;
};

#endif
