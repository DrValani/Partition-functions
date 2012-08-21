#ifndef AddLayers2_CC
#define AddLayers2_CC

#include <vector>
#include "Globals.cc"
#include "Storage.cc"
#include "refPolynomial.cc"
#include "commonMethods.cc"

int join_spins_to_base(const Storage &vlist, const Storage &GrowthVector,
    vector<refPolynomial> &vrefpol) {
	refPolynomial refpol2;
	addEnergy4pol(vlist, GrowthVector, refpol2, MAXPOWER);
	vrefpol.push_back(refpol2);
	return vrefpol.size() - 1;
}

int create_new_vlist(int current_column, const Storage &GrowthVector,
    const int &columnsAdded, vector<vector<Storage> >  &vlist, const int vlistIndex,
    vector<refPolynomial> &vrefpol, vector<vector<Storage> >  &vpol_vlist) {

	if (current_column < ySize - 1) {
		Storage souter(true);
		for (int config = 0; config < ARRSIZE; config++) {
			int vector_index = create_new_vlist(current_column + 1, GrowthVector,
			    columnsAdded, vlist, vlist[columnsAdded][vlistIndex][config], vrefpol,
			    vpol_vlist);
			souter.push_back(config, vector_index);

		}
		if (current_column + 1 == ySize - 1) {
			vpol_vlist[columnsAdded + 1].push_back(souter);
			return vpol_vlist[columnsAdded + 1].size() - 1;
		} else {
			vlist[columnsAdded + 1].push_back(souter);
			return vlist[columnsAdded + 1].size() - 1;
		}
	} else {
		return join_spins_to_base(vpol_vlist[columnsAdded][vlistIndex],
		    GrowthVector, vrefpol);

	}
}

void CalculatePartitionVector(int current_column,
    const G_TABLE1 &current_config, G_TABLE1 &previous_config,
    vector<Storage>   &GrowthMatrix, vector<vector<Storage> >  &vlist, vector<vector<Polynomial> > &vpol,
    const int &vpolnew, vector<vector<Storage> >  &vpol_vlist, vector<vector<refPolynomial> >  &vrefpol,
    Polynomial &p) {
	if (current_column < ySize - 1) {

		previous_config[current_column] = current_config[current_column];
		int configuration = current_config[current_column];  // The fixed configuration of the extra column of spin
		//Calculate new polynomial
		const int columnsAdded = current_column;
		//vpol[current_column + 2].clear();
		vlist[current_column + 1].clear();
		vpol_vlist[current_column + 1].clear();
		vrefpol[current_column].clear();

		//findmaxpower(GrowthMatrix[configuration]);
		int vlist_index = vlist[current_column].size() - 1;

		create_new_vlist(current_column, GrowthMatrix[configuration], columnsAdded,
		    vlist, vlist_index, vrefpol[current_column], vpol_vlist);

		if (vrefpol[current_column].size() < MAX_VPOL_SIZE) {
			calculate_inner_vpol(current_column, vpol, vrefpol, p);
			vrefpol[current_column].clear();
		}
		previous_config[current_column + 1] = -1;
	} else {
		refPolynomial refpol;
		int vlistIndex = vpol_vlist[current_column].size() - 1;
		int configuration = current_config[current_column];

		addEnergy4pol(vpol_vlist[current_column][vlistIndex],
		    GrowthMatrix[configuration], refpol, MAXPOWER);

		int Hamiltonian = calculate_outer_hamiltonian(GrowthMatrix, current_config);
		vpol[0][vpolnew].clear();
		refpol.sum(vpol[0][vpolnew], vpol[current_column + 1], MAXPOWER,
		    Hamiltonian, p);
	}
}

/*
 *This method will scroll through the ylist -ie the list of configurations
 */

void ScrollThroughAllConfigurations(const G_TABLE2 &ylist, vector<vector<Storage> >  &vlist,
    vector<vector<Polynomial> > &vpol, vector<vector<Storage> >  &vpol_vlist, vector<vector<refPolynomial> >  &vrefpol,
    Polynomial &p) {

	vector<Storage>   GrowthMatrix;
	Create_energy_list(GrowthMatrix);
	Storage s;
	Create_energy_list(s);
	GrowthMatrix.push_back(s);
	vector<int> previous_config(ySize + 1, -1);

	for (int config = 0; config < ylist.size(); config++) {  //This will scroll through all configurations

		int current_column = 0;
		while (current_column < ySize) {
			while (ylist[config][current_column] == previous_config[current_column])
				current_column++;
			CalculatePartitionVector(current_column, ylist[config], previous_config,
			    GrowthMatrix, vlist, vpol, ylist[config][ySize], vpol_vlist, vrefpol,
			    p);
			current_column++;
		}
	}
}

#endif

