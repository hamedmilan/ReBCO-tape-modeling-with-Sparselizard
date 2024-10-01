//
//
// This headers sums the strain and stress within a loop of loadings
// it uses custom function!
//
//
std::vector<densemat> sum_strain_stress(std::vector<densemat> vals, std::vector<field> infields,
                                    elementselector& elsel,        std::vector<double>& gps,
                                    expression* meshdeform)
{
    //
    // rows & cols of the dense matrices
    int nrows = vals[0].countrows();
    int ncols = vals[0].countcolumns();

    // declaring pointers for strain_e from input
    double* e_xx = vals[0].getvalues();
    double* e_yy = vals[1].getvalues();
    double* e_zz = vals[2].getvalues();
    double* e_yz = vals[3].getvalues();
    double* e_xz = vals[4].getvalues();
    double* e_xy = vals[5].getvalues();

    // declaring pointers for stress_e from input
    double* s_xx = vals[6].getvalues();
    double* s_yy = vals[7].getvalues();
    double* s_zz = vals[8].getvalues();
    double* s_yz = vals[9].getvalues();
    double* s_xz = vals[10].getvalues();
    double* s_xy = vals[11].getvalues();

    // declaring pointers for strain_tot from input (exp of strain_tot)
    double* e_tot_xx = vals[12].getvalues();
    double* e_tot_yy = vals[13].getvalues();
    double* e_tot_zz = vals[14].getvalues();
    double* e_tot_yz = vals[15].getvalues();
    double* e_tot_xz = vals[16].getvalues();
    double* e_tot_xy = vals[17].getvalues();

    // declaring pointers for stress_corr from input (exp of stress_corr)
    double* sig_corr_xx = vals[18].getvalues();
    double* sig_corr_yy = vals[19].getvalues();
    double* sig_corr_zz = vals[20].getvalues();
    double* sig_corr_yz = vals[21].getvalues();
    double* sig_corr_xz = vals[22].getvalues();
    double* sig_corr_xy = vals[23].getvalues();

    // going through all nodes
    for (int i = 0; i < nrows*ncols; i++) {


        //std::cout << i << '\n';
        // updating strain tensor
        e_tot_xx[i] = e_tot_xx[i] + e_xx[i];
        e_tot_yy[i] = e_tot_yy[i] + e_yy[i];
        e_tot_zz[i] = e_tot_zz[i] + e_zz[i];
        //
        e_tot_yz[i] = e_tot_yz[i] + e_yz[i];
        e_tot_xz[i] = e_tot_xz[i] + e_xz[i];
        e_tot_xy[i] = e_tot_xy[i] + e_xy[i];
        /*
        e_tot_yz[i] = e_tot_yz[i] + 0.5* e_yz[i];
        e_tot_xz[i] = e_tot_xz[i] + 0.5* e_xz[i];
        e_tot_xy[i] = e_tot_xy[i] + 0.5* e_xy[i];
        */

        // updating stress tensor
        sig_corr_xx[i] = sig_corr_xx[i] + s_xx[i];
        sig_corr_yy[i] = sig_corr_yy[i] + s_yy[i];
        sig_corr_zz[i] = sig_corr_zz[i] + s_zz[i];
        //
        sig_corr_yz[i] = sig_corr_yz[i] + s_yz[i];
        sig_corr_xz[i] = sig_corr_xz[i] + s_xz[i];
        sig_corr_xy[i] = sig_corr_xy[i] + s_xy[i];

    }

    // writing new strain_tot to its field
    infields[0].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[12]);
    infields[1].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[13]);
    infields[2].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[14]);
    infields[3].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[15]);
    infields[4].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[16]);
    infields[5].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[17]);
    // writing new stress_corr to its field
    infields[6].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[18]);
    infields[7].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[19]);
    infields[8].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[20]);
    infields[9].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[21]);
    infields[10].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[22]);
    infields[11].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[23]);



    // returning some dummy value!
    return {vals[12]};
}
