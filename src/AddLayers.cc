#ifndef AddLayers_CC
#define AddLayers_CC

#include <iostream>
#include "latticeSize.cc"
#include "Globals.cc"
#include "Polynomial.cc"
#include "Storage.cc"
#include "commonMethods.cc"
#include "fileMethods.cc"
#include "myTimer.cc"
#include "build2d_9.cc"
#include "build2d_10.cc"
#include "AddLayers2.cc"
using namespace std;

int join_spins_to_base(
    const G_TABLE1 &transformOperation,  // Probably gonna need to transform vlist before application
    Storage vlist, const Storage &GrowthVector, vector<refPolynomial> &vrefpol,
    G_TABLE2 &hst2, const int &vlistIndex) {
	int gp = 0;
	for (int i = 0; i < transformOperation.size(); i++) {
		gp = GROUPTABLE[gp][transformOperation[i]];
	}
	if (gp != 0) vlist = TRANSFORMGROUP[gp].applyGroup(vlist);
	int pin = hst2[vlistIndex][gp];
	if (pin == -1) {
		refPolynomial refpol2;
		addEnergy4pol(vlist, GrowthVector, refpol2, MAXPOWER);
		hst2[vlistIndex][gp] = search_add_list(vrefpol, refpol2, hst2[vlistIndex]);
	}
	return hst2[vlistIndex][gp];
}

/**
 * This method will be used to find the vector associated with @param-config.
 */
int find_vlist_vector(const Storage &vlist, const int &config,
    const G_TABLE1 &vlist_transformGroup, G_TABLE1 &sub_vlist_transformGroup,
    G_TABLE1 &transformOperations) {
	int new_config = config;
	//for (int i = transformOperations.size() - 1; i >= 0; i--){
	for (int i = 0; i < transformOperations.size(); i++) {
		new_config = TRANSFORMGROUP[transformOperations[i]][new_config];
	}

	if (vlist[new_config] == -1) {
		//Then  we got to find and apply its symmetry
		int operationIndex = -1;
		int sym_config = inverse_find_symmetrical_equivalent(new_config,
		    vlist_transformGroup, vlist, operationIndex);
		if (operationIndex == -1 || sym_config == -1) {
			std::cout << "on no" << std::endl;
			std::cout << "config-> " << config << " sym_config->" << sym_config
			    << std::endl;
			std::cout << vlist << std::endl;
			exit(0);
		}
		Calculate_subgroup(sub_vlist_transformGroup, sym_config,
		    vlist_transformGroup);
		transformOperations.push_back(
		    TRANSFORMGROUP[vlist_transformGroup[operationIndex]].getId());
		return vlist[sym_config];
	} else {
		Calculate_subgroup(sub_vlist_transformGroup, new_config,
		    vlist_transformGroup);
		return vlist[new_config];
	}
}

int create_new_vlist(int current_column, const G_TABLE1 transformGroup,
    const G_TABLE1 vlist_transformGroup, const Storage &GrowthVector,
    const int &columnsAdded, vector<vector<Storage> >  &vlist, const int vlistIndex,
    G_TABLE1 transformOperations, vector<refPolynomial> &vrefpol,
    vector<vector<Storage> >  &vpol_vlist, G_TABLE2 &hst2) {

	if (current_column < ySize - 1) {
		Storage souter(true);
		for (int config = 0; config < ARRSIZE; config++) {
			int sym_equivalent_config = find_symmetrical_equivalent(config,
			    transformGroup, souter);
			if (sym_equivalent_config == -1) {

				G_TABLE1 sub_transformGroup, sub_vlist_transformGroup,
				    sub_transformOperations = transformOperations;
				Calculate_subgroup(sub_transformGroup, config, transformGroup);
				//expands out old vlist
				int vlistIndex_new = find_vlist_vector(vlist[columnsAdded][vlistIndex],
				    config, vlist_transformGroup, sub_vlist_transformGroup,
				    sub_transformOperations);
				int vector_index = create_new_vlist(current_column + 1,
				    sub_transformGroup, sub_vlist_transformGroup, GrowthVector,
				    columnsAdded, vlist, vlistIndex_new, sub_transformOperations,
				    vrefpol, vpol_vlist, hst2);
				souter.push_back(config, vector_index);
			} else {
				if (current_column + 1 == ySize - 1) souter.push_back(config,
				    souter[sym_equivalent_config]);
				else souter.push_back(config, -1);
			}

		}
		if (current_column + 1 == ySize - 1) return search_add_list_exact(
		    vpol_vlist[columnsAdded + 1], souter);
		else {
			vlist[columnsAdded + 1].push_back(souter);
			return vlist[columnsAdded + 1].size() - 1;
		}
	} else {
		return join_spins_to_base(transformOperations,
		    vpol_vlist[columnsAdded][vlistIndex], GrowthVector, vrefpol, hst2,
		    vlistIndex);

	}
}

void create_hash_table(G_TABLE2 &hst2, int vpol_vlist_size) {
	hst2.clear();
	G_TABLE1 temp(TRANSFORMGROUP.size(), -1);
	for (int i = 0; i < vpol_vlist_size; i++) {
		hst2.push_back(temp);
	}
}

void CalculatePartitionVector(int current_column,
    const G_TABLE1 &current_config, G_TABLE1 &previous_config,
    vector<Storage>   &GrowthMatrix, G_TABLE2 &transformGroups,
    G_TABLE2 &vlist_transformGroup, vector<vector<Storage> >  &vlist, vector<vector<Polynomial> > &vpol,
    const int &vpolnew, vector<vector<Storage> >  &vpol_vlist, G_TABLE3 & hst3,
    vector<vector<refPolynomial> >  &vrefpol, Polynomial &p) {
	if (current_column < ySize - 1) {

		previous_config[current_column] = current_config[current_column];
		int configuration = current_config[current_column];  // The fixed configuration of the extra column of spin
		//Calculate new polynomial
		transformGroups[current_column + 1].clear();
		Calculate_subgroup(transformGroups[current_column + 1], configuration,
		    transformGroups[current_column]);
		const int columnsAdded = current_column;
		vlist[current_column + 1].clear();
		vpol_vlist[current_column + 1].clear();
		vrefpol[current_column].clear();

		G_TABLE1 transformOperations;

		int vlist_index = vlist[current_column].size() - 1;
		create_hash_table(hst3[current_column], vpol_vlist[current_column].size());

		create_new_vlist(current_column, transformGroups[current_column + 1],
		    vlist_transformGroup[current_column], GrowthMatrix[configuration],
		    columnsAdded, vlist, vlist_index, transformOperations,
		    vrefpol[current_column], vpol_vlist, hst3[current_column]);

		if (vrefpol[current_column].size() < MAX_VPOL_SIZE) {
			calculate_inner_vpol(current_column, vpol, vrefpol, p);
			vrefpol[current_column].clear();
			SMALLEST_TRANSFERMATRIX = current_column;
		} else {
			cout << "vrefpol[" << current_column << "].size() = "
			    << vrefpol[current_column].size() << endl;
			//SMALLEST_TRANSFERMATRIX=0;
		}
		vlist_transformGroup[current_column + 1].clear();
		Calculate_subgroup(vlist_transformGroup[current_column + 1], configuration,
		    vlist_transformGroup[current_column]);

		previous_config[current_column + 1] = -1;
	} else {
		refPolynomial refpol;
		int vlistIndex = vpol_vlist[current_column].size() - 1;
		int configuration = current_config[current_column];

		addEnergy4pol(vpol_vlist[current_column][vlistIndex],
		    GrowthMatrix[configuration], refpol, MAXPOWER);

		int Hamiltonian = calculate_outer_hamiltonian(GrowthMatrix, current_config);
		if (vpol[1].size() < MAX_VPOL_SIZE) {
			vpol[0][vpolnew].clear();
			refpol.sum(vpol[0][vpolnew], vpol[current_column + 1], MAXPOWER,
			    Hamiltonian, p);
		} else {
			vpol[0][0].clear();
			refpol.sum(vpol[0][0], vpol[current_column + 1], MAXPOWER, Hamiltonian,
			    p);
			pol_to_vpol_file(VPOL_FILENAME, vpol[0][0]);
		}
	}
}
void CalculatePartitionVectorOnly(int current_column,
    const G_TABLE1 &current_config, G_TABLE1 &previous_config,
    vector<Storage>   &GrowthMatrix, G_TABLE2 &transformGroups,
    G_TABLE2 &vlist_transformGroup, vector<vector<Storage> >  &vlist, vector<vector<Polynomial> > &vpol,
    const int &vpolnew, vector<vector<Storage> >  &vpol_vlist, G_TABLE3 & hst3,
    vector<vector<refPolynomial> >  &vrefpol, Polynomial &p, const int &start_column) {
	if (current_column < ySize - 1) {

		//previous_config[current_column] = current_config[current_column];
		int configuration = current_config[current_column];  // The fixed configuration of the extra column of spin
		//Calculate new polynomial
		transformGroups[current_column + 1].clear();
		Calculate_subgroup(transformGroups[current_column + 1], configuration,
		    transformGroups[current_column]);
		const int columnsAdded = current_column;
		vlist[current_column + 1].clear();
		vpol_vlist[current_column + 1].clear();
		vrefpol[current_column].clear();

		G_TABLE1 transformOperations;

		int vlist_index = vlist[current_column].size() - 1;
		create_hash_table(hst3[current_column], vpol_vlist[current_column].size());

		create_new_vlist(current_column, transformGroups[current_column + 1],
		    vlist_transformGroup[current_column], GrowthMatrix[configuration],
		    columnsAdded, vlist, vlist_index, transformOperations,
		    vrefpol[current_column], vpol_vlist, hst3[current_column]);

		vlist_transformGroup[current_column + 1].clear();
		Calculate_subgroup(vlist_transformGroup[current_column + 1], configuration,
		    vlist_transformGroup[current_column]);

		previous_config[current_column + 1] = -1;
	} else {
		vector<vector<refPolynomial> >  vrefpolnew(ySize + 1);
		int vlistIndex = vpol_vlist[current_column].size() - 1;
		int configuration = current_config[current_column];
		refPolynomial refpol;
		addEnergy4pol(vpol_vlist[current_column][vlistIndex],
		    GrowthMatrix[configuration], refpol, MAXPOWER);
		vrefpol[ySize - 1].push_back(refpol);
		// sum up refpols
		int mxpow = 1;
		int cc = 0;
		int starter = 0;
		if (SMALLEST_TRANSFERMATRIX < start_column) cc = SMALLEST_TRANSFERMATRIX;
		lookagain: while (cc < start_column && vrefpol[cc].size() == 0) {
			cc++;
		}
		int cc1 = cc + 1;
		while (cc1 < start_column && vrefpol[cc1].size() == 0) {
			cc = cc1;
			goto lookagain;
		}
		bool SWAP = true;
		if (cc < start_column && vrefpol[cc].size() > 0) {
			vrefpolnew[cc].swap(vrefpol[cc]);
			starter = cc;
			while (cc < ySize - 1) {
				for (int i = 0; i < vrefpol[cc + 1].size(); i++) {
					refpol.clear();
					refpol.sum(vrefpol[cc + 1][i], vrefpolnew[cc], MAXPOWER * mxpow,
					    MAXPOWER);
					vrefpolnew[cc + 1].push_back(refpol);
				}
				if (SWAP) vrefpolnew[cc].swap(vrefpol[cc]);
				SWAP = false;
				mxpow++;
				cc++;
			}
		} else {
			vrefpolnew[start_column].swap(vrefpol[start_column]);
			starter = start_column;
			cc = start_column;
			while (cc < ySize - 1) {
				for (int i = 0; i < vrefpol[cc + 1].size(); i++) {
					refpol.clear();
					refpol.sum(vrefpol[cc + 1][i], vrefpolnew[cc], MAXPOWER * mxpow,
					    MAXPOWER);
					vrefpolnew[cc + 1].push_back(refpol);
				}
				if (SWAP) vrefpolnew[start_column].swap(vrefpol[start_column]);
				SWAP = false;
				mxpow++;
				cc++;
			}
		}

		int Hamiltonian = calculate_outer_hamiltonian(GrowthMatrix, current_config);
		if (vpol[1].size() < MAX_VPOL_SIZE) {
			vpol[0][vpolnew].clear();
			vrefpolnew[ySize - 1][0].sum(vpol[0][vpolnew], vpol[starter + 1],
			    MAXPOWER * mxpow, Hamiltonian, p);
		} else {
			vpol[0][0].clear();
			vrefpolnew[ySize - 1][0].sum(vpol[0][0], vpol[starter + 1],
			    MAXPOWER * mxpow, Hamiltonian, p);
			pol_to_vpol_file(VPOL_FILENAME, vpol[0][0]);
		}
		vrefpol[ySize - 1].clear();
	}
}

void ScrollThroughAllConfigurations(const G_TABLE2 &ylist,
    G_TABLE1 &transformGroup, vector<vector<Storage> >  &vlist, vector<vector<Polynomial> > &vpol,
    vector<vector<Storage> >  &vpol_vlist, vector<vector<refPolynomial> >  &vrefpol, Polynomial &p) {

	G_TABLE2 transformGroup2(ySize + 1), vlist_transformGroup2(ySize + 1);
	transformGroup2[0] = transformGroup;
	vlist_transformGroup2[0] = transformGroup;
	vector<Storage>   GrowthMatrix;
	Create_energy_list(GrowthMatrix);
	Storage s;
	Create_energy_list(s);
	GrowthMatrix.push_back(s);
	vector<int> previous_config(ySize + 1, -1);
	G_TABLE3 hst3(ySize + 1);
	time_t startTime_inner;
	time(&startTime_inner);
	for (int config = 0; config < ylist.size(); config++) {  //This will scroll through all configurations

		int current_column = 0;
		while (current_column < ySize) {
			while (ylist[config][current_column] == previous_config[current_column])
				current_column++;
			if (config < ylist.size() - 1
			    && ylist[config][current_column]
			        != ylist[config + 1][current_column]) {
				if (current_column < ySize - 2) cout << flush << current_column << " is on its own " << ylist[config] << "/" << ylist.size()
				    << "\r";
				int columns_added = current_column;
				while (current_column < ySize) {
					CalculatePartitionVectorOnly(current_column, ylist[config],
					    previous_config, GrowthMatrix, transformGroup2,
					    vlist_transformGroup2, vlist, vpol, ylist[config][ySize],
					    vpol_vlist, hst3, vrefpol, p, columns_added);
					current_column++;
				}
			} else {
				CalculatePartitionVector(current_column, ylist[config], previous_config,
				    GrowthMatrix, transformGroup2, vlist_transformGroup2, vlist, vpol,
				    ylist[config][ySize], vpol_vlist, hst3, vrefpol, p);
				cout << flush << "Progress in current layer " << ylist[config] << config
				    << "/" << ylist.size() << "\r";

				if (current_column == 1) {
					cout << endl;
					calculateETA(startTime_inner, config + 1, ylist.size());
				}
				current_column++;
			}
		}
	}
}

void buildAdd(G_TABLE2 &ylist, G_TABLE1 &transformGroup, vector<vector<Storage> >  &vlist,
    vector<vector<Polynomial> > &vpol, vector<vector<Storage> >  &vpol_vlist, const refPolynomial & total) {

	std::cout << "Now running build add " << std::endl;

	vector<vector<refPolynomial> >  vrefpol(ySize + 1);

	Polynomial p;
	p.clear();
	total.sum(p, vpol[1]);
	std::cout << "Printing p\n";
	cout, p;
	cout << endl;
	print2File(p, 1);
	final_check(p, 1);
	if (vpol[1].size() < MAX_VPOL_SIZE) vpol[0] = vpol[1];
	else {
		truncateFile (VPOL_FILENAME);
		vpol[0].push_back(p);
	}
	int totalCalculations = zSize - 1;
	int calculated = 0;

	for (int z = 1; z < zSize; z++) {
		cout << "Now adding layer " << z + 1 << endl;
		ScrollThroughAllConfigurations(ylist, transformGroup, vlist, vpol,
		    vpol_vlist, vrefpol, p);

		if (vpol[1].size() < MAX_VPOL_SIZE) vpol[0].swap(vpol[1]);
		else {
			file_to_vpol(vpol[1], ylist, VPOL_FILENAME);
			if (z != zSize - 1) truncateFile (VPOL_FILENAME);
		}

		p.clear();
		total.sum(p, vpol[1]);
		std::cout << "Printing p\n";
		cout, p;
		cout << endl;
		print2File(p, z + 1);
		final_check(p, z + 1);
		std::cout << "just made " << xSize << "x" << ySize << "x" << z + 1
		    << std::endl;
		std::cout << zSize - z - 1 << " blocks left to add" << std::endl;
		calculated++;
		cout << "OVERALL ";
		calculateETA(startTime, calculated, totalCalculations);
	}
	for (int i = 0; i < ySize + 1; i++) {
		cout << "vpol[" << i << "].size = " << vpol[i].size() << endl;
	}
}

void buildAdd(G_TABLE2 &ylist_bc3, G_TABLE1 &transformGroup_bc3,
    vector<vector<Storage> >  &vlist, vector<vector<Polynomial> > &vpol, vector<vector<Storage> >  & vpol_vlist) {

	vector<vector<refPolynomial> >  vrefpol(ySize + 1);
	Polynomial p;
	std::cout << "Now running build add " << std::endl;
	vector<Storage>   GrowthMatrix;
	Create_energy_matrix(GrowthMatrix);
	G_TABLE2 ylist;
	int totalCalculations = zSize * ylist_bc3.size();
	int calculated = 0;
	vector<Polynomial> vpol_final(zSize);
	for (int i = 0; i < zSize; i++)
		vpol_final[i].clear();
	for (int bc3_config = 0; bc3_config < ylist_bc3.size(); bc3_config++) {
		G_TABLE1 transformGroup = transformGroup_bc3;
		for (int i = 0; i < ySize; i++) {
			G_TABLE1 subTransformGroup;
			Calculate_subgroup(subTransformGroup, ylist_bc3[bc3_config][i],
			    transformGroup);
			transformGroup = subTransformGroup;
		}
		//clearLists(vpol);
		clearLists(vpol_vlist);
		clearLists(vlist);
		ylist.clear();
		build_2dLattice(vlist, transformGroup, ylist, vpol[1], vpol_vlist[0],
		    ylist_bc3[bc3_config], GrowthMatrix);  //Make the 2d lattice. (8)
		int z_vpol_ref = find_ylist_reference(ylist, ylist_bc3[bc3_config]);
		vpol[0] = vpol[1];
		vpol_final[0].multiPlus(vpol[1][z_vpol_ref], ylist_bc3[bc3_config][ySize]);
		cout << "-------------REMAINING " << ylist_bc3.size() - bc3_config;
		cout << " Layers added " << 1 << endl;
		calculated++;
		calculateETA(startTime, calculated, totalCalculations);
		if (transformGroup.size() > 0) {
			for (int z = 1; z < zSize; z++) {
				if (z == zSize - 1) {
					ylist.clear();
					ylist.push_back(ylist_bc3[bc3_config]);
					ylist[0][ySize] = z_vpol_ref;
				}
				ScrollThroughAllConfigurations(ylist, transformGroup, vlist, vpol,
				    vpol_vlist, vrefpol, p);
				vpol_final[z].multiPlus(vpol[0][z_vpol_ref],
				    ylist_bc3[bc3_config][ySize]);
				vpol[0].swap(vpol[1]);
				//ostringstream out;
				cout << "-------------REMAINING " << ylist_bc3.size() - bc3_config;
				cout << " Layers added " << z + 1 << endl;
				//cout << out;
				calculated++;
				calculateETA(startTime, calculated, totalCalculations);
			}
		} else {
			for (int z = 1; z < zSize; z++) {
				if (z == zSize - 1) {
					ylist.clear();
					ylist.push_back(ylist_bc3[bc3_config]);
					ylist[0][ySize] = z_vpol_ref;
				}
				ScrollThroughAllConfigurations(ylist, vlist, vpol, vpol_vlist, vrefpol,
				    p);
				vpol_final[z].multiPlus(vpol[0][z_vpol_ref],
				    ylist_bc3[bc3_config][ySize]);
				vpol[0].swap(vpol[1]);
				cout << "-------------REMAINING " << ylist_bc3.size() - bc3_config;
				cout << " Layers added " << z + 1 << endl;
				calculated++;
				calculateETA(startTime, calculated, totalCalculations);
			}
		}
	}
	for (int z = 0; z < zSize; z++) {
		std::cout << "Printing " << z + 1 << "\n";
		print2File(vpol_final[z], z + 1);
		std::cout, vpol_final[z];
		std::cout << std::endl;
		final_check(vpol_final[z], z + 1);
	}
}

#endif

