#include "xdock_header.h"

int main(int argc, char** argv) {
	int idur;

	auto told = std::chrono::high_resolution_clock::now();
	int ires;
	InputParameters *ip;
	XData *xd;
	ip = new InputParameters();
	ires = ip->readInputParameters(argc, argv);
	if (ires != 0) {
		delete ip; ip = nullptr;
		return -1;
	}

	xd = new XData(ip);

	ires = xd->startProcessing();
	
	auto tnow = std::chrono::high_resolution_clock::now();
	idur = std::chrono::duration_cast<std::chrono::milliseconds>(tnow - told).count();
	if (xd->flog != nullptr) {
		fprintf(xd->flog, "\nTime spent (in s): %i.%03i\n", idur / 1000, idur % 1000);
	}
	printf("\nTime spent (in s): %i.%03i\n", idur / 1000, idur % 1000);
		
	delete xd; xd = nullptr;
	delete ip; ip = nullptr;

	return ires;
}
