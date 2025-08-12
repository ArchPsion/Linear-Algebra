#ifndef __LINEAR_ALGEBRA_HPP__
#define __LINEAR_ALGEBRA_HPP__

#include "Using.hpp"
#include <array>
#include <iostream>

template <typename Type>
Type Unpad(Type&& str)
{
	const auto pos = str.find_first_of('.');
	
	if (pos != std::string::npos)
	{
		while (str.back() == '0')
			str.pop_back();
		
		if (str.back() == '.')
			str.pop_back();		
	}
	
	return str;
}

class HexAlgebra
{
	public:
	
		template <typename Type>
		static Type Conjugate(const Type&);
};

template <>
f32 HexAlgebra::Conjugate(const f32& f)
{
	return f;
}

template <>
f64 HexAlgebra::Conjugate(const f64& f)
{
	return f;
}

template <typename Type>
Type HexAlgebra::Conjugate(const Type& c)
{
	return c.conj();
}

template <typename, u32, u32>
class HexMatrix;

template <typename Type, u32 Number>
class HexVector
{
	private:
	
		std::array<Type, Number> 			coordinates;
	
	public:
	
		template <typename... OtherType>
		static HexVector<Type, Number> 			Make(OtherType...);
		template <typename OtherType>
		static HexVector<Type, Number> 			Make(OtherType&&);
		
		static HexVector<Type, Number> 			Zero(void);
		
		template <typename... OtherType>
		HexVector(OtherType...);
		template <typename OtherType>
		HexVector(OtherType&&);
		
		Type						scalar(const HexVector<Type, Number>&) const;
		void						show(void) const;
		
		template <u32 NumberOfColumns>
		HexVector<Type, NumberOfColumns>		operator*(const HexMatrix<Type, Number, NumberOfColumns>&) const;
		template <u32 NumberOfColumns>
		HexMatrix<Type, Number, NumberOfColumns>	operator*(const HexVector<Type, NumberOfColumns>&) const;
		Type&						operator()(u32);
		const Type&					operator()(u32) const;
		
		template <typename, u32>
		friend class HexVector;
		
		template <typename, u32, u32>
		friend class HexMatrix;
};

template <typename Type, u32 Number>
template <typename... OtherType>
HexVector<Type, Number>::HexVector(OtherType... args) : coordinates({ args... })
{
	static_assert(Number != 0u);
}

template <typename Type, u32 Number>
template <typename OtherType>
HexVector<Type, Number>::HexVector(OtherType&& container)
{
	static_assert(Number != 0u);
	
	const auto end = container.cend();
	auto it = container.cbegin();
	
	for (auto& c : coordinates)
	{
		if (it == end)
			break;
		
		c = *it;
		++it;
	}
}

template <typename Type, u32 Number>
template <typename... OtherType>
HexVector<Type, Number> HexVector<Type, Number>::Make(OtherType... args)
{
	return HexVector<Type, Number>(args...);
}

template <typename Type, u32 Number>
template <typename OtherType>
HexVector<Type, Number> HexVector<Type, Number>::Make(OtherType&& container)
{
	return HexVector<Type, Number>(container);
}

template <typename Type, u32 Number>
HexVector<Type, Number> HexVector<Type, Number>::Zero(void)
{
	return HexVector<Type, Number>();
}

template <typename Type, u32 Number>
template <u32 NumberOfColumns>
HexVector<Type, NumberOfColumns> HexVector<Type, Number>::operator*(const HexMatrix<Type, Number, NumberOfColumns>& m) const
{
	HexVector<Type, NumberOfColumns> result;
	auto it = result.coordinates.begin();
	
	for (auto i = 0u; i < NumberOfColumns; ++i)
	{
		*it = m.scalarColumn(i, *this);
		++it;
	}
	
	return result;
}

template <typename Type, u32 Number>
template <u32 NumberOfColumns>
HexMatrix<Type, Number, NumberOfColumns> HexVector<Type, Number>::operator*(const HexVector<Type, NumberOfColumns>& v) const
{
	HexMatrix<Type, Number, NumberOfColumns> result;
	auto it = result.values.begin();
	
	for (const auto& s : HexVector::coordinates)
	{
		for (const auto& t : v.coordinates)
		{
			*it = s*t;
			++it;
		}
	}
	
	return result;
}

template <typename Type, u32 Number>
Type HexVector<Type, Number>::scalar(const HexVector<Type, Number>& v) const
{
	auto sum = static_cast<Type>(0.f);
	auto it = v.coordinates.cbegin();
	
	for (const auto& s : HexVector::coordinates)
	{
		sum += HexAlgebra::Conjugate(*it)*s;
		++it;
	}
	
	return sum;
}

template <typename Type, u32 Number>
void HexVector<Type, Number>::show(void) const
{
	std::cout << "\t[,1]";
	
	auto it = HexVector::coordinates.cbegin();
	std::cout << std::endl;
	
	for (auto i = 0u; i < Number; ++i)
	{
		std::cout << '[' << i+1u << ",]\t" << ::Unpad(std::to_string(*it)) << std::endl;
		++it;
	}
}

template <typename Type, u32 Number>
Type& HexVector<Type, Number>::operator()(u32 index)
{
	return HexVector::coordinates[index];
}

template <typename Type, u32 Number>
const Type& HexVector<Type, Number>::operator()(u32 index) const
{
	return HexVector::coordinates[index];
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
class HexMatrix
{
	private:
	
		std::array<Type, NumberOfRows*NumberOfColumns>		values;
		
		Type							detPivot(void);
		HexMatrix<Type, NumberOfRows, NumberOfColumns>		inversePivot(void);
	
	public:
	
		template <typename... OtherType>
		static HexMatrix<Type, NumberOfRows, NumberOfColumns> 	Make(OtherType...);
		template <typename OtherType>
		static HexMatrix<Type, NumberOfRows, NumberOfColumns> 	Make(OtherType&&);
		
		static HexMatrix<Type, NumberOfRows, NumberOfColumns>	Identity(void);
		static HexMatrix<Type, NumberOfRows, NumberOfColumns> 	Zero(void);
		
		template <typename OtherType>
		HexMatrix(OtherType&&);
		template <typename... OtherType>
		HexMatrix(OtherType...);
		
		void							addColumn(u32, Type, u32, u32 = 0u);
		void							addRow(u32, Type, u32, u32 = 0u);
		void							clear(void);
		Type							det(void) const;
		HexMatrix<Type, NumberOfRows, NumberOfColumns>		inverse(void) const;
		void							multiplyColumn(u32, Type, u32 = 0u);
		void							multiplyRow(u32, Type, u32 = 0u);
		HexMatrix<Type, NumberOfRows, NumberOfColumns>		power(u32) const;
		Type							scalar(const HexMatrix<Type, NumberOfRows, NumberOfColumns>&) const;
		Type							scalarColumn(u32, const HexVector<Type, NumberOfRows>&) const;
		Type							scalarRow(u32, const HexVector<Type, NumberOfColumns>&) const;
		void							show(void) const;
		void							swapColumns(u32, u32, u32 = 0u);
		void							swapRows(u32, u32, u32 = 0u);
		Type							trace(void) const;
		HexMatrix<Type, NumberOfColumns, NumberOfRows>		transpose(void) const;
		
		template <u32 Number>
		HexMatrix<Type, NumberOfRows, Number>			operator*(const HexMatrix<Type, NumberOfColumns, Number>&) const;
		HexVector<Type, NumberOfRows>				operator*(const HexVector<Type, NumberOfColumns>&) const;
		void							operator*=(const HexMatrix<Type, NumberOfRows, NumberOfColumns>&);
		Type&							operator()(u32, u32);
		const Type&						operator()(u32, u32) const;
		
		template <typename, u32>
		friend class HexVector;
		
		template <typename, u32, u32>
		friend class HexMatrix;
};

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
template <typename... OtherType>
HexMatrix<Type, NumberOfRows, NumberOfColumns>::HexMatrix(OtherType... args) : values({ args... })
{
	static_assert(NumberOfRows != 0u);
	static_assert(NumberOfColumns != 0u);
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
template <typename OtherType>
HexMatrix<Type, NumberOfRows, NumberOfColumns>::HexMatrix(OtherType&& container)
{
	static_assert(NumberOfRows != 0u);
	static_assert(NumberOfColumns != 0u);
	
	const auto end = container.cend();
	auto it = container.cbegin();
	
	for (auto& v : values)
	{
		if (it == end)
			break;
		
		v = *it;
		++it;
	}
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
HexMatrix<Type, NumberOfRows, NumberOfColumns> HexMatrix<Type, NumberOfRows, NumberOfColumns>::Identity(void)
{
	static_assert(NumberOfRows == NumberOfColumns);
	
	HexMatrix<Type, NumberOfRows, NumberOfColumns> result;
	auto it = result.values.begin();
	
	for (auto i = 0u; i < NumberOfColumns; ++i)
	{
		*it = static_cast<Type>(1.f);
		it += (NumberOfColumns + 1u);
	}
	
	return result;
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
template <typename... OtherType>
HexMatrix<Type, NumberOfRows, NumberOfColumns> HexMatrix<Type, NumberOfRows, NumberOfColumns>::Make(OtherType... args)
{
	return HexMatrix<Type, NumberOfRows, NumberOfColumns>(args...);
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
template <typename OtherType>
HexMatrix<Type, NumberOfRows, NumberOfColumns> HexMatrix<Type, NumberOfRows, NumberOfColumns>::Make(OtherType&& container)
{
	return HexMatrix<Type, NumberOfRows, NumberOfColumns>(container);
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
HexMatrix<Type, NumberOfRows, NumberOfColumns> HexMatrix<Type, NumberOfRows, NumberOfColumns>::Zero(void)
{
	return HexMatrix<Type, NumberOfRows, NumberOfColumns>();
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
template <u32 Number>
HexMatrix<Type, NumberOfRows, Number> HexMatrix<Type, NumberOfRows, NumberOfColumns>::operator*(const HexMatrix<Type, NumberOfColumns, Number>& m) const
{
	HexMatrix<Type, NumberOfRows, Number> result;
	auto it = result.values.begin();
	
	for (auto i = 0u; i < NumberOfRows; ++i)
	{
		for (auto j = 0u; j < Number; ++j)
		{
			auto it1 = HexMatrix::values.cbegin() + i*NumberOfColumns;
			auto it2 = m.values.cbegin() + j;
			
			for (auto k = 0u; k < NumberOfColumns; ++k)
			{
				*it += (*it1)*(*it2);
				++it1;
				it2 += Number;
			}
			
			++it;
		}
	}
	
	return result;
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
HexVector<Type, NumberOfRows> HexMatrix<Type, NumberOfRows, NumberOfColumns>::operator*(const HexVector<Type, NumberOfColumns>& v) const
{
	HexVector<Type, NumberOfRows> result;
	
	for (auto i = 0u; i < NumberOfRows; ++i)
		result(i) = HexMatrix::scalarRow(i, v);
	
	return result;
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
void HexMatrix<Type, NumberOfRows, NumberOfColumns>::operator*=(const HexMatrix<Type, NumberOfRows, NumberOfColumns>& m)
{
	static_assert(NumberOfRows == NumberOfColumns);
	HexMatrix::values = HexMatrix::operator*(m).values;
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
Type& HexMatrix<Type, NumberOfRows, NumberOfColumns>::operator()(u32 row, u32 column)
{
	return HexMatrix::values[row*NumberOfColumns + column];
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
const Type& HexMatrix<Type, NumberOfRows, NumberOfColumns>::operator()(u32 row, u32 column) const
{
	return HexMatrix::values[row*NumberOfColumns + column];
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
void HexMatrix<Type, NumberOfRows, NumberOfColumns>::addColumn(u32 c1, Type factor, u32 c2, u32 r)
{
	auto it1 = HexMatrix::values.begin() + r*NumberOfColumns + c1;
	auto it2 = HexMatrix::values.cbegin() + r*NumberOfColumns + c2;
	
	while (r < NumberOfRows)
	{
		*it1 += factor*(*it2);
		it1 += NumberOfColumns;
		it2 += NumberOfColumns;
		++r;
	}
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
void HexMatrix<Type, NumberOfRows, NumberOfColumns>::addRow(u32 r1, Type factor, u32 r2, u32 c)
{
	auto it1 = HexMatrix::values.begin() + r1*NumberOfColumns + c;
	auto it2 = HexMatrix::values.cbegin() + r2*NumberOfColumns + c;
	
	while (c < NumberOfColumns)
	{
		*it1 += factor*(*it2);
		++it1;
		++it2;
		++c;
	}
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
void HexMatrix<Type, NumberOfRows, NumberOfColumns>::clear(void)
{
	HexMatrix::values.fill(static_cast<Type>(0.f));
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
Type HexMatrix<Type, NumberOfRows, NumberOfColumns>::det(void) const
{
	static_assert(NumberOfRows == NumberOfColumns);
	
	if constexpr (NumberOfColumns < 3u)
		return (NumberOfColumns < 2u ? HexMatrix::values.back() : HexMatrix::values[0u]*HexMatrix::values[3u] - HexMatrix::values[1u]*HexMatrix::values[2u]);
	else
	{
		auto copy = *this;
		return copy.detPivot();
	}
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
Type HexMatrix<Type, NumberOfRows, NumberOfColumns>::detPivot(void)
{
	auto result = static_cast<Type>(1.f);
	auto swaps = 0u;
	
	for (auto i = 0u; i < NumberOfRows - 1u; ++i)
	{
		auto pivotRow = i;
		
		while (pivotRow < NumberOfRows and HexMatrix::operator()(pivotRow, i) == static_cast<Type>(0.f))
			++pivotRow;
		
		if (pivotRow > i)
		{
			if (pivotRow < NumberOfRows)
			{
				++swaps;
				HexMatrix::swapRows(i, pivotRow, i);
			}
			else
				return static_cast<Type>(0.f);
		}
		
		const auto div = HexMatrix::operator()(i, i);
		result *= div;
		
		for (auto j = i + 1u; j < NumberOfRows; ++j)
		{
			const auto factor = -HexMatrix::operator()(j, i)/div;
			
			if (static_cast<Type>(0.f) != factor)
				HexMatrix::addRow(j, factor, i, i + 1u);
		}
	}
	
	return (swaps % 2u != 0u ? -result : result)*HexMatrix::values.back();
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
HexMatrix<Type, NumberOfRows, NumberOfColumns> HexMatrix<Type, NumberOfRows, NumberOfColumns>::inverse(void) const
{
	static_assert(NumberOfRows == NumberOfColumns);
	
	if constexpr (NumberOfColumns < 3u)
	{
		const auto det = HexMatrix::det();
		
		if (static_cast<Type>(0.f) != det)
			return (NumberOfColumns < 2u ? HexMatrix::Make(static_cast<Type>(1.f/det)) : HexMatrix::Make(HexMatrix::values[3u]/det, -HexMatrix::values[1u]/det, -HexMatrix::values[2u]/det, HexMatrix::values[0u]));
		
		return HexMatrix::Zero();
	}
	else
	{
		auto copy = *this;
		return copy.inversePivot();
	}
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
HexMatrix<Type, NumberOfRows, NumberOfColumns> HexMatrix<Type, NumberOfRows, NumberOfColumns>::inversePivot(void)
{
	auto result = HexMatrix<Type, NumberOfRows, NumberOfColumns>::Identity();
	
	for (auto i = 0u; i < NumberOfRows; ++i)
	{
		auto pivotRow = i;
		
		while (pivotRow < NumberOfRows and HexMatrix::operator()(pivotRow, i) == static_cast<Type>(0.f))
			++pivotRow;
		
		if (pivotRow > i)
		{
			if (pivotRow < NumberOfRows)
			{
				HexMatrix::swapRows(i, pivotRow, i);
				result.swapRows(i, pivotRow);
			}
			else
			{
				result.values.fill(static_cast<Type>(0.f));
				break;
			}
		}
		
		const auto div = HexMatrix::operator()(i, i);
		HexMatrix::multiplyRow(i, static_cast<Type>(1.f)/div, i);
		result.multiplyRow(i, static_cast<Type>(1.f)/div);
		
		for (auto j = 0u; j < NumberOfRows; ++j)
		{
			if (j == i)
				continue;
			
			const auto factor = -HexMatrix::operator()(j, i);
			
			if (static_cast<Type>(0.f) != factor)
			{
				HexMatrix::addRow(j, factor, i);
				result.addRow(j, factor, i);
			}
		}
	}
	
	return result;
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
void HexMatrix<Type, NumberOfRows, NumberOfColumns>::multiplyColumn(u32 c, Type factor, u32 r)
{
	auto it = HexMatrix::values.begin() + NumberOfColumns*r + c;
	
	while (r < NumberOfRows)
	{
		*it *= factor;
		it += NumberOfColumns;
		++r;
	}
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
void HexMatrix<Type, NumberOfRows, NumberOfColumns>::multiplyRow(u32 r, Type factor, u32 c)
{
	auto it = HexMatrix::values.begin() + NumberOfColumns*r + c;
	
	while (c < NumberOfColumns)
	{
		*it *= factor;
		++it;
		++c;
	}
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
HexMatrix<Type, NumberOfRows, NumberOfColumns> HexMatrix<Type, NumberOfRows, NumberOfColumns>::power(u32 p) const
{
	switch (p)
	{
		case 0u:
			return HexMatrix<Type, NumberOfRows, NumberOfColumns>::Identity();
		
		case 1u:
			return *this;
		
		default:
			break;
	}
	
	if (p & 1u)
		return HexMatrix::operator*(HexMatrix::operator*(*this).power(p >> 1u));
	
	return HexMatrix::operator*(*this).power(p >> 1u);
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
Type HexMatrix<Type, NumberOfRows, NumberOfColumns>::scalar(const HexMatrix<Type, NumberOfRows, NumberOfColumns>& m) const
{
	auto sum = static_cast<Type>(0.f);
	auto it = m.values.cbegin();
	
	for (const auto& s : HexMatrix::values)
	{
		sum += HexAlgebra::Conjugate(*it)*s;
		++it;
	}
	
	return sum;
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
Type HexMatrix<Type, NumberOfRows, NumberOfColumns>::scalarColumn(u32 column, const HexVector<Type, NumberOfRows>& v) const
{
	auto sum = static_cast<Type>(0.f);
	auto it = HexMatrix::values.cbegin() + column;
	
	for (const auto& s : v.coordinates)
	{
		sum += (*it)*s;
		it += NumberOfColumns;
	}
	
	return sum;
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
Type HexMatrix<Type, NumberOfRows, NumberOfColumns>::scalarRow(u32 row, const HexVector<Type, NumberOfColumns>& v) const
{
	auto sum = static_cast<Type>(0.f);
	auto it = HexMatrix::values.cbegin() + row*NumberOfColumns;
	
	for (const auto& s : v.coordinates)
	{
		sum += (*it)*s;
		++it;
	}
	
	return sum;
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
void HexMatrix<Type, NumberOfRows, NumberOfColumns>::show(void) const
{
	std::cout << "\t[.1]";
	
	for (auto i = 1u; i < NumberOfColumns; ++i)
		std::cout << "\t\t[." << i+1u << ']';
	
	auto it = HexMatrix::values.cbegin();
	std::cout << std::endl;
	
	for (auto i = 0u; i < NumberOfRows; ++i)
	{
		auto str = ::Unpad(std::to_string(*it));
		++it;
		
		std::cout << '[' << i+1u << ".]\t" << str;
		
		for (auto j = 1u; j < NumberOfColumns; ++j)
		{
			if (str.size() < 8u)
				std::cout << '\t';
			
			str = ::Unpad(std::to_string(*it));
			++it;
			
			std::cout << '\t' << str;
		}
		
		std::cout << std::endl;
	}
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
void HexMatrix<Type, NumberOfRows, NumberOfColumns>::swapColumns(u32 c1, u32 c2, u32 r)
{
	auto it1 = HexMatrix::values.begin() + r*NumberOfColumns + c1;
	auto it2 = HexMatrix::values.begin() + r*NumberOfColumns + c2;
	
	while (r < NumberOfRows)
	{
		std::swap(*it1, *it2);
		it1 += NumberOfColumns;
		it2 += NumberOfColumns;
		++r;
	}
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
void HexMatrix<Type, NumberOfRows, NumberOfColumns>::swapRows(u32 r1, u32 r2, u32 c)
{
	auto it1 = HexMatrix::values.begin() + r1*NumberOfColumns + c;
	auto it2 = HexMatrix::values.begin() + r2*NumberOfColumns + c;
	
	while (c < NumberOfColumns)
	{
		std::swap(*it1, *it2);
		++it1;
		++it2;
		++c;
	}
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
Type HexMatrix<Type, NumberOfRows, NumberOfColumns>::trace(void) const
{
	static_assert(NumberOfRows == NumberOfColumns);
	auto it = HexMatrix::values.cbegin();
	auto sum = *it;
	
	for (auto k = 0u; k < NumberOfColumns - 1u; ++k)
	{
		it += (NumberOfColumns + 1u);
		sum += *it;
	}
	
	return sum;
}

template <typename Type, u32 NumberOfRows, u32 NumberOfColumns>
HexMatrix<Type, NumberOfColumns, NumberOfRows> HexMatrix<Type, NumberOfRows, NumberOfColumns>::transpose(void) const
{
	HexMatrix<Type, NumberOfColumns, NumberOfRows> result;
	auto it = HexMatrix::values.cbegin();
	auto count = 0u;
	
	for (auto& v : result.values)
	{
		v = *it;
		++count;
		it = (count % NumberOfRows != 0u ? it + NumberOfColumns : HexMatrix::values.cbegin() + count/NumberOfRows);
	}
	
	return result;
}

#endif
