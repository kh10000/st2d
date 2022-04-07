#include <math.h>
#include <stdio.h>
#include "calculate_chebyshev.h"

extern "C" int calculate_chebyshev(double** cell, double** cart, double** scale,
                                   int* atom_i, double* atom_weight, int natoms, int* cal_atoms, int cal_num,
                                   int* params_i, double* params_d, int nsyms,
                                   double** symf) {
    // cell: cell info of structure
    // cart: cartesian coordinates of atoms
    // scale: fractional coordinates of atoms
    // atom_i: atom type index (start with 1)
    // atom_weight: weight assigned for each atom type
    // params_i: order of expansion
    //           [N_rad, N_ang, N_rad(weightd), N_ang(weighted)]
    // params_d: cutoff distances
    //           [Rc for rdf, Rc for adf, Rc for weighted-rdf, Rc for weighted-adf]
    // natoms: # of atoms
    // nsyms: descriptor dimension
    // symf: descriptor vector ([# of atoms, # of symfuncs])

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, cutoff, dradtmp, rRij, rRik, rRjk;
    double plane_d[3], total_shift[3], precal[5], tmpd[9], dangtmp[3], lcoeff[9];
    double vecij[3], vecik[3], vecjk[3], deljk[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];
    double max_rc_ang = 0.0;
    double xij, xijk;
    // int nsf[5+1];

    int rad_dim   = params_i[0]+1;
    int ang_dim   = params_i[1]+1;
    int rad_dim_w = params_i[2]+1;
    int ang_dim_w = params_i[3]+1;
 
    double *p_rad   = new double[rad_dim];
    double *p_ang   = new double[ang_dim];
    double *p_rad_w = new double[rad_dim_w];
    double *p_ang_w = new double[ang_dim_w];
    p_rad[0] = 1.0; p_ang[0] = 1.0; p_rad_w[0] = 1.0; p_ang_w[0] = 1.0;

    // double *u_rad   = new double[params_i[0]];
    // double *u_ang   = new double[params_i[1]];
    // double *u_rad_w = new double[params_i[2]];
    // double *u_ang_w = new double[params_i[3]];
    // u_rad[0] = 1.0; u_ang[0] = 1.0; u_rad_w[0] = 1.0; u_ang_w[0] = 1.0;

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    for (int s=1; s < 4; s+=2) {
        max_rc_ang = max_rc_ang > params_d[s] ? max_rc_ang : params_d[s];
    }

    cutoff = 0.0;
    for (int s=0; s < 4; ++s) {
        if (cutoff < params_d[s])
            cutoff = params_d[s];
    }

    total_bins = 1;

    // calculate the inverse matrix of cell and the distance between cell plane
    cross[0][0] = cell[1][1]*cell[2][2] - cell[1][2]*cell[2][1];
    cross[0][1] = cell[1][2]*cell[2][0] - cell[1][0]*cell[2][2];
    cross[0][2] = cell[1][0]*cell[2][1] - cell[1][1]*cell[2][0];
    cross[1][0] = cell[2][1]*cell[0][2] - cell[2][2]*cell[0][1];
    cross[1][1] = cell[2][2]*cell[0][0] - cell[2][0]*cell[0][2];
    cross[1][2] = cell[2][0]*cell[0][1] - cell[2][1]*cell[0][0];
    cross[2][0] = cell[0][1]*cell[1][2] - cell[0][2]*cell[1][1];
    cross[2][1] = cell[0][2]*cell[1][0] - cell[0][0]*cell[1][2];
    cross[2][2] = cell[0][0]*cell[1][1] - cell[0][1]*cell[1][0];

    vol = cross[0][0]*cell[0][0] + cross[0][1]*cell[0][1] + cross[0][2]*cell[0][2];
    inv[0][0] = cross[0][0]/vol;
    inv[0][1] = cross[1][0]/vol;
    inv[0][2] = cross[2][0]/vol;
    inv[1][0] = cross[0][1]/vol;
    inv[1][1] = cross[1][1]/vol;
    inv[1][2] = cross[2][1]/vol;
    inv[2][0] = cross[0][2]/vol;
    inv[2][1] = cross[1][2]/vol;
    inv[2][2] = cross[2][2]/vol;

    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
        }
        bin_i[i][3] = bin_i[i][0] + nbins[0]*bin_i[i][1] + nbins[0]*nbins[1]*bin_i[i][2];
        atoms_bin[bin_i[i][3]]++;
    }

    max_atoms_bin = 0;
    for (int i=0; i < total_bins; ++i) {
        if (atoms_bin[i] > max_atoms_bin)
            max_atoms_bin = atoms_bin[i];
    }

    delete[] atoms_bin;

    // # of bins in each direction
    neigh_check_bins = 1;
    for (int i=0; i < 3; ++i) {
        bin_range[i] = ceil(cutoff * nbins[i] / plane_d[i]);
        neigh_check_bins *= 2*bin_range[i];
    }

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 5 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        nneigh = 0;

        for (int j=0; j < 3; ++j) {
            max_bin[j] = bin_i[i][j] + bin_range[j];
            min_bin[j] = bin_i[i][j] - bin_range[j];
        }

        for (int dx=min_bin[0]; dx < max_bin[0]+1; ++dx) {
            for (int dy=min_bin[1]; dy < max_bin[1]+1; ++dy) {
                for (int dz=min_bin[2]; dz < max_bin[2]+1; ++dz) {
                    pbc_bin[0] = (dx%nbins[0] + nbins[0]) % nbins[0];
                    pbc_bin[1] = (dy%nbins[1] + nbins[1]) % nbins[1];
                    pbc_bin[2] = (dz%nbins[2] + nbins[2]) % nbins[2];
                    cell_shift[0] = (dx-pbc_bin[0]) / nbins[0];
                    cell_shift[1] = (dy-pbc_bin[1]) / nbins[1];
                    cell_shift[2] = (dz-pbc_bin[2]) / nbins[2];

                    bin_num = pbc_bin[0] + nbins[0]*pbc_bin[1] + nbins[0]*nbins[1]*pbc_bin[2];

                    for (int j=0; j < natoms; ++j) {
                        if (bin_i[j][3] != bin_num)
                            continue;

                        if (!(cell_shift[0] || cell_shift[1] || cell_shift[2]) && (i == j))
                            continue;

                        for (int a=0; a < 3; ++a) {
                            total_shift[a] = cell_shift[0]*cell[0][a] + cell_shift[1]*cell[1][a] + cell_shift[2]*cell[2][a]
                                             + cart[j][a] - cart[i][a];
                        }

                        tmp = sqrt(total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2]);

                        if (tmp < cutoff) {
                            for (int a=0; a < 3; ++a)
                                nei_list_d[nneigh*5 + a] = total_shift[a];
                            nei_list_d[nneigh*5 + 3] = tmp;
                            nei_list_d[nneigh*5 + 4] = atom_weight[j];
                            nei_list_i[nneigh*2]    = atom_i[j];
                            nei_list_i[nneigh*2 + 1] = j;
                            nneigh++;
                        }
                    }
                }
            }
        }

        for (int j=0; j < nneigh; ++j) {
            // calculate expansion cofficients of RDF and weighted-RDF
            // RDF
            rRij = nei_list_d[j*5 + 3];

            xij = 2 * rRij/params_d[0] - 1;
            precal[0] = cutf(rRij/params_d[0]);
            // precal[1] = dcutf(rRij, params_d[0]);
            
            if (rad_dim > 0){
                p_rad[1] = xij;
                // u_rad[1] = 2*xij;
            }

            for (int l=2; l < rad_dim; ++l){
                p_rad[l] = 2*xij*p_rad[l-1] - p_rad[l-2];
                // u_rad[l] = 2*xij*u_rad[l-1] - u_rad[l-2];
            }

            for (int s=0; s < rad_dim; ++s) {
                symf[ii][s] += p_rad[s]*precal[0];
            }

            // weighted rdf
            xij = 2*rRij/params_d[2] - 1;
            precal[0] = cutf(rRij/params_d[2]);
            // precal[1] = dcutf(rRij, params_d[2]);

            if (rad_dim_w > 0){
                p_rad_w[1] = xij;
                // u_rad_w[1] = 2*xij;
            }

            for (int l=2; l < rad_dim_w; ++l){
                p_rad_w[l] = 2*xij*p_rad_w[l-1] - p_rad_w[l-2];
                // u_rad_w[l] = 2*xij*u_rad_w[l-1] - u_rad_w[l-2];
            }

            for (int s=0; s < rad_dim_w; ++s){
                symf[ii][s+rad_dim+ang_dim] += p_rad_w[s]*precal[0]*nei_list_d[j*5+4];
            }

            if (rRij > max_rc_ang) continue;
            for (int k=j+1; k < nneigh; ++k) {
                // ADF
                rRik = nei_list_d[k*5 + 3];
                if (rRik > max_rc_ang) continue;

                deljk[0] = nei_list_d[k*5]     - nei_list_d[j*5];
                deljk[1] = nei_list_d[k*5 + 1] - nei_list_d[j*5 + 1];
                deljk[2] = nei_list_d[k*5 + 2] - nei_list_d[j*5 + 2];
                rRjk = sqrt(deljk[0]*deljk[0] + deljk[1]*deljk[1] + deljk[2]*deljk[2]);

                if (rRjk < 0.0001) continue;

                precal[0] = cutf(rRij/params_d[1]);
                precal[1] = cutf(rRik/params_d[1]);
                precal[2] = cutf(rRij/params_d[3]);
                precal[3] = cutf(rRik/params_d[3]);
                precal[4]  = (rRij*rRij + rRik*rRik - rRjk*rRjk)/2/rRij/rRik;

                xijk  = precal[4];
                if (ang_dim > 0) {
                    p_ang[1] = xijk;
                    // u_ang[1] = 2*xijk;
                }
                
                for (int l=2; l < ang_dim; ++l){
                    p_ang[l] = 2*xijk*p_ang[l-1] - p_ang[l-2];
                    // u_ang[l] = 2*xijk*u_ang[l-1] - u_ang[l-2];
                }

                for (int s=0; s < ang_dim; ++s) {
                    symf[ii][s+rad_dim] += p_ang[s]*precal[0]*precal[1];                    
                }

                // weighted ADF
                if (ang_dim_w > 0){
                    p_ang_w[1] = xijk;
                    // u_ang_w[1] = 2*xijk;
                }

                for (int l=2; l < ang_dim_w; ++l){
                    p_ang_w[l] = 2*xijk*p_ang_w[l-1] - p_ang[l-2];
                    // u_ang_w[l] = 2*xijk*u_ang_w[l-1] - u_ang[l-2];
                }

                for (int s=0; s < ang_dim_w; ++s){
                    symf[ii][s+rad_dim+ang_dim+rad_dim_w] += \
                            p_ang_w[s]*precal[2]*precal[3]*nei_list_d[5*j+4]*nei_list_d[5*k+4];
                }
            }
        }

        delete[] nei_list_d;
        delete[] nei_list_i;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;

    delete[] p_rad;
    delete[] p_ang;
    delete[] p_rad_w;
    delete[] p_ang_w;

    return 0;
}

void PyInit_libsymf(void) { } // for windows
