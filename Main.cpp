#include "LinearAlgebra.hpp"
#include <vector>

#include <csignal>
#include <execinfo.h>
#include <unistd.h>

void Handler(i32)
{
	void* tab[20u];
	
	const auto size = backtrace(tab, 20u);
	backtrace_symbols_fd(tab, size, STDERR_FILENO);
	std::exit(1);
}

i32 main(void)
{
	struct sigaction act;
	sigemptyset(&act.sa_mask);
	act.sa_handler = ::Handler;
	act.sa_flags = 0;
	sigaction(SIGSEGV, &act, 0);
	
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
