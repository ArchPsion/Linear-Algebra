// Custom Libraries
#include "LinearAlgebra.hpp"

int main(void)
{
	const auto m1 = HexMatrix<f32, 2u, 2u>(1.f, 1.f, 0.f, 1.f);
	const auto m2 = HexMatrix<f32, 2u, 3u>(2.f, 3.f, 4.f, 5.f, 6.f, 7.f);
	const auto m3 = HexMatrix<f32, 3u, 3u>(1.f, 0.f, 0.f, 2.f, 1.f, 1.f, 0.f, 6.f, 0.f);
	
	const auto v1 = HexVector<f32, 2u>(20.f, 40.f);
	const auto v2 = HexVector<f32, 3u>(5.f, 10.f, 15.f);
	
	const auto v3 = v1*m1;
	const auto v4 = v1*m2;
	
	m3.show();
	m3.inverse().show();
	return 0;
}
