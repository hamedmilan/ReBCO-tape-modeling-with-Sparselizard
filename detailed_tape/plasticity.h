//
//
// This headers : plasticity CF
// it uses custom function!
//
//
std::vector<densemat> plasticity(std::vector<densemat> vals,
                                    std::vector<field> infields,
                                    elementselector& elsel,
                                    std::vector<double>& gps,
                                    expression* meshdeform);
// relative stress
std::vector<double> relative_stress(std::vector<double> S, std::vector<double> alpha);
// yield function
double yield_function(std::vector<double> ksi, double k);
// deviatoric sigma
std::vector<double> dev(std::vector<double> t);
// n_hat
std::vector<double> n_hat(std::vector<double> st);
// delta Lambda
double delta_lambda(double F, double mu, double H);
// L2norm
double L2norm(std::vector<double> t);
// *****************************************************************************
// *****************************************************************************
std::vector<densemat> plasticity(std::vector<densemat> vals, std::vector<field> infields,
                                    elementselector& elsel,        std::vector<double>& gps,
                                    expression* meshdeform)
{
    // plasticity function
    // ****************************
    // ****************************

    // general plasticity idea
    // checks if the integration point is in plastic regime
    // if it is
    // (1) corrects the stress           : stress_corr
    // (2) calculates the plastic strain : strain_p
    // (3) evolves the back stress       : alpha
    // (4) evolves the sigma_y           : kappa/sigma_y
    // (5) computes the correction term for Consistent Tangent Operator (CTO)
    // !!! by CTO we mean the correction term of it
    // if not
    // (1) adds incremental strain to total strain
    // (2) adds incremental stress to total stress

    // !!! this function is triggered through plasticity_trigger formulation
    // !!! so, it will be called 3 times!
    // !!! here we check that, if plasticity is called before
    // if no, we do the plasticity
    // if yes, we skip!
    /*
    call_id++;
    if (call_id > 1) {

        if (call_id == 3) {
            // if 3, we reset the call_id
            call_id = 0;
        }

        // in case of 2, or 3
        return {vals[15], vals[16], vals[17],
                vals[18], vals[19], vals[20]};
    }
    */

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

    // declaring pointers for mu, H & beta from input
    double* mu = vals[12].getvalues();
    double* H = vals[13].getvalues();
    double* beta = vals[14].getvalues();

    // declaring pointers for strain_tot from input (exp of strain_tot)
    double* e_tot_xx = vals[15].getvalues();
    double* e_tot_yy = vals[16].getvalues();
    double* e_tot_zz = vals[17].getvalues();
    double* e_tot_yz = vals[18].getvalues();
    double* e_tot_xz = vals[19].getvalues();
    double* e_tot_xy = vals[20].getvalues();

    // declaring pointers for stress_corr from input (exp of stress_corr)
    double* sig_corr_xx = vals[21].getvalues();
    double* sig_corr_yy = vals[22].getvalues();
    double* sig_corr_zz = vals[23].getvalues();
    double* sig_corr_yz = vals[24].getvalues();
    double* sig_corr_xz = vals[25].getvalues();
    double* sig_corr_xy = vals[26].getvalues();

    // declaring pointers for strain_p from input (exp of strain_p)
    double* e_p_xx = vals[27].getvalues();
    double* e_p_yy = vals[28].getvalues();
    double* e_p_zz = vals[29].getvalues();
    double* e_p_yz = vals[30].getvalues();
    double* e_p_xz = vals[31].getvalues();
    double* e_p_xy = vals[32].getvalues();

    // declaring pointers for alpha from input (exp of alpha)
    double* a_xx = vals[33].getvalues();
    double* a_yy = vals[34].getvalues();
    double* a_zz = vals[35].getvalues();
    double* a_yz = vals[36].getvalues();
    double* a_xz = vals[37].getvalues();
    double* a_xy = vals[38].getvalues();

    // declaring pointer for sigma_y from input (exp of sigma_y)
    double* kpp = vals[39].getvalues();

    // declaring pointers for C_CTO from input (exp of CTO)
    double* cto_11 = vals[40].getvalues();
    double* cto_12 = vals[41].getvalues();
    double* cto_13 = vals[42].getvalues();
    double* cto_14 = vals[43].getvalues();
    double* cto_15 = vals[44].getvalues();
    double* cto_16 = vals[45].getvalues();

    double* cto_21 = vals[46].getvalues();
    double* cto_22 = vals[47].getvalues();
    double* cto_23 = vals[48].getvalues();
    double* cto_24 = vals[49].getvalues();
    double* cto_25 = vals[50].getvalues();
    double* cto_26 = vals[51].getvalues();

    double* cto_31 = vals[52].getvalues();
    double* cto_32 = vals[53].getvalues();
    double* cto_33 = vals[54].getvalues();
    double* cto_34 = vals[55].getvalues();
    double* cto_35 = vals[56].getvalues();
    double* cto_36 = vals[57].getvalues();

    double* cto_41 = vals[58].getvalues();
    double* cto_42 = vals[59].getvalues();
    double* cto_43 = vals[60].getvalues();
    double* cto_44 = vals[61].getvalues();
    double* cto_45 = vals[62].getvalues();
    double* cto_46 = vals[63].getvalues();

    double* cto_51 = vals[64].getvalues();
    double* cto_52 = vals[65].getvalues();
    double* cto_53 = vals[66].getvalues();
    double* cto_54 = vals[67].getvalues();
    double* cto_55 = vals[68].getvalues();
    double* cto_56 = vals[69].getvalues();

    double* cto_61 = vals[70].getvalues();
    double* cto_62 = vals[71].getvalues();
    double* cto_63 = vals[72].getvalues();
    double* cto_64 = vals[73].getvalues();
    double* cto_65 = vals[74].getvalues();
    double* cto_66 = vals[75].getvalues();

    // declaring pointers for x,y & z
    double* x = vals[76].getvalues();
    double* y = vals[77].getvalues();
    double* z = vals[78].getvalues();
    double* p_update_flag = vals[79].getvalues();


    // node-ish values

    // temporary values for easy handling
    std::vector<double> sig_(6);    // stress tensor
    std::vector<double> ep_(6);     // plastic strain tensor
    std::vector<double> alpha_(6);  // back stress tensor
    double kappa_;                  // sigma_y

    // plasticity loop params.
    std::vector<double> ksi(6);     // relative stress : dev_sig - alpha
    std::vector<double> nHat(6);    // n_hat : the unit vector of ksi_trial
    double F;                       // yield_function
    double dLambda;                 // delta Lambda : the absolute change of strain_p


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

        // temporary stress
        sig_ = {sig_corr_xx[i],
                sig_corr_yy[i],
                sig_corr_zz[i],
                sig_corr_yz[i],
                sig_corr_xz[i],
                sig_corr_xy[i]};

        // temporary plastic strain
        ep_ = {e_p_xx[i], e_p_yy[i], e_p_zz[i], e_p_yz[i], e_p_xz[i], e_p_xy[i]};

        // temporary alpha
        alpha_ = {a_xx[i], a_yy[i], a_zz[i], a_yz[i], a_xz[i], a_xy[i]};

        // temporary kappa
        kappa_ = kpp[i];

        // making CTO (temporary CTO)
        /*
        densemat cto_i(6,6, {cto_11[i], cto_12[i], cto_13[i], cto_14[i], cto_15[i], cto_16[i],
                             cto_21[i], cto_22[i], cto_23[i], cto_24[i], cto_25[i], cto_26[i],
                             cto_31[i], cto_32[i], cto_33[i], cto_34[i], cto_35[i], cto_36[i],
                             cto_41[i], cto_42[i], cto_43[i], cto_44[i], cto_45[i], cto_46[i],
                             cto_51[i], cto_52[i], cto_53[i], cto_54[i], cto_55[i], cto_56[i],
                             cto_61[i], cto_62[i], cto_63[i], cto_64[i], cto_65[i], cto_66[i] });
        */

        // initializing the CTO to zero (we are not using the previous CTO values!)
        densemat cto_i(6,6, 0.0);
        //cto_i.print();

        //keeeping kappa0 for whole loop
        //kappa0_ = kpp[i];
        //kappa_ = kappa0_;

        // relative stress
        ksi = relative_stress(dev(sig_), alpha_);

        // yield function
        F = yield_function(ksi, kappa_);

        // n_hat
        nHat = n_hat(ksi);

        // initializing delta_lambda
        dLambda = 0.0;

        // CTO stuff!
        densemat n_outer_n;         // n_hat (outer product) n_hat :: (n_hat*n_hat)
        densemat proj_dev;          // deviatoric projection matrix
        densemat dev_nn(6,6, 0.0);  // subtraction of n_hat*n_hat from dev.proj.

        // plasticity flag
        int p_flag = 0;

        // if plasticity, compute:
        // CTO stuff
        if (F-1e-7>0) {

            // changing the plasticity flag!
            p_flag = 1;

            // re-arranging n_hat in 3*3 mat form
            double nhat_mat[3][3] =
            {
            {nHat[0], nHat[5], nHat[4]},
            {nHat[5], nHat[1], nHat[3]},
            {nHat[4], nHat[3], nHat[2]}
            };

            // declaring n_hat*n_hat in the form of 4th order tensor (3*3*3*3)
            double N_N[3][3][3][3];

            // computing n_hat*n_hat
            for (int ii = 0; ii < 3; ii++) {
                for (int jj = 0; jj < 3; jj++) {
                    //

                    for (int kk = 0; kk < 3; kk++) {
                        for (int ll = 0; ll < 3; ll++) {
                            //
                            N_N[ii][kk][ll][jj] = nhat_mat[ii][kk] * nhat_mat[ll][jj];
                        }

                    }


                }
            }

            // re-arranging n_hat*n_hat in 6*6 mat form
            n_outer_n = densemat(6,6, {N_N[0][0][0][0], N_N[0][0][1][1], N_N[0][0][2][2], N_N[0][0][1][2], N_N[0][0][0][2], N_N[0][0][0][1],
                                       N_N[0][0][1][1], N_N[1][1][1][1], N_N[1][1][2][2], N_N[1][1][1][2], N_N[1][1][0][2], N_N[1][1][0][1],
                                       N_N[0][0][2][2], N_N[1][1][2][2], N_N[2][2][2][2], N_N[2][2][1][2], N_N[2][2][0][2], N_N[2][2][0][1],
                                       N_N[0][0][1][2], N_N[1][1][1][2], N_N[2][2][1][2], N_N[1][2][1][2], N_N[1][2][0][2], N_N[1][2][0][1],
                                       N_N[0][0][0][2], N_N[1][1][0][2], N_N[2][2][0][2], N_N[1][2][0][2], N_N[0][2][0][2], N_N[0][2][0][1],
                                       N_N[0][0][0][1], N_N[1][1][0][1], N_N[2][2][0][1], N_N[1][2][0][1], N_N[0][2][0][1], N_N[0][1][0][1]});

            // deviatoric projection in 6*6 mat form
            /*
            proj_dev = densemat(6,6, {1.0-1.0/3.0,    -1.0/3.0,    -1.0/3.0, 0.0, 0.0, 0.0,
                                         -1.0/3.0, 1.0-1.0/3.0,    -1.0/3.0, 0.0, 0.0, 0.0,
                                         -1.0/3.0,    -1.0/3.0, 1.0-1.0/3.0, 0.0, 0.0, 0.0,
                                              0.0,         0.0,         0.0, 1.0, 0.0, 0.0,
                                              0.0,         0.0,         0.0, 0.0, 1.0, 0.0,
                                              0.0,         0.0,         0.0, 0.0, 0.0, 1.0,});*/
            proj_dev = densemat(6,6, {1.0-1.0/3.0,    -1.0/3.0,    -1.0/3.0, 0.0, 0.0, 0.0,
                                         -1.0/3.0, 1.0-1.0/3.0,    -1.0/3.0, 0.0, 0.0, 0.0,
                                         -1.0/3.0,    -1.0/3.0, 1.0-1.0/3.0, 0.0, 0.0, 0.0,
                                              0.0,         0.0,         0.0, 0.5, 0.0, 0.0,
                                              0.0,         0.0,         0.0, 0.0, 0.5, 0.0,
                                              0.0,         0.0,         0.0, 0.0, 0.0, 0.5,});

            // making dev.proj. - n_hat*n_hat
            dev_nn = proj_dev.copy();
            dev_nn.subtract(n_outer_n);

            // making the first term of CTO subtraction
            densemat mat1 = n_outer_n.getproduct((4.0*mu[i]*mu[i])/(2.0*mu[i] + H[i]));
            cto_i.add(mat1);
            //mat1.print();
            //std::cout << '\n';

        }

        double D_lambda = 0.0;

        // plasticity loop
        int ploopcount = 0;
        while (F-1e-7>0) {
            //std::cout << "It goes to plastic!" << "\n\n";
            //abort();

            // if plasticity reaches a max number!
            if (ploopcount>1000) {
                //std::cout << "!!! F = " << F << '\n';
                break;
            }

            // computing dLambda and updating stuff
            dLambda = delta_lambda(F, mu[i], H[i]);

            // updating plastic strain
            ep_[0] = ep_[0] + dLambda* nHat[0];
            ep_[1] = ep_[1] + dLambda* nHat[1];
            ep_[2] = ep_[2] + dLambda* nHat[2];
            //
            ep_[3] = ep_[3] + 2.0* dLambda* nHat[3];
            ep_[4] = ep_[4] + 2.0* dLambda* nHat[4];
            ep_[5] = ep_[5] + 2.0* dLambda* nHat[5];

            // updating stress
            sig_[0] = sig_[0] - 2.0*mu[i]*dLambda* nHat[0];
            sig_[1] = sig_[1] - 2.0*mu[i]*dLambda* nHat[1];
            sig_[2] = sig_[2] - 2.0*mu[i]*dLambda* nHat[2];
            //
            sig_[3] = sig_[3] - 2.0*mu[i]*dLambda* nHat[3];
            sig_[4] = sig_[4] - 2.0*mu[i]*dLambda* nHat[4];
            sig_[5] = sig_[5] - 2.0*mu[i]*dLambda* nHat[5];

            // updating kappa
            kappa_ = kappa_ + beta[i] * H[i] * dLambda;

            D_lambda = D_lambda + dLambda;

            // updating CTO
            densemat mat2 = dev_nn.getproduct((4.0*dLambda*mu[i]*mu[i])/L2norm(ksi));
            cto_i.add(mat2);

            //mat2.print();
            //abort();

            // updating ksi & F based on new stress
            ksi = relative_stress(dev(sig_), alpha_);
            F = yield_function(ksi, kappa_);

            // updating n_hat (needed?!)
            nHat = n_hat(ksi);

            ploopcount++;
        }


        // writing temporary values to pointer (if plasticity happened!)
        if (p_flag == 1) {

            //densemat mat2 = dev_nn.getproduct((4.0*D_lambda*mu[i]*mu[i])/L2norm(ksi));
            //cto_i.add(mat2);

            // temp_stress >> stress tensor
            sig_corr_xx[i] = sig_[0];
            sig_corr_yy[i] = sig_[1];
            sig_corr_zz[i] = sig_[2];
            //
            sig_corr_yz[i] = sig_[3];
            sig_corr_xz[i] = sig_[4];
            sig_corr_xy[i] = sig_[5];

            // temp_ep >> plastic strain
            e_p_xx[i] = ep_[0];
            e_p_yy[i] = ep_[1];
            e_p_zz[i] = ep_[2];
            //
            e_p_yz[i] = ep_[3];
            e_p_xz[i] = ep_[4];
            e_p_xy[i] = ep_[5];

            // temp_kappa >> kappa
            kpp[i] = kappa_;


            // temp_CTO >> CTO
            cto_11[i] = cto_i.getvalue(0, 0);
            cto_12[i] = cto_i.getvalue(0, 1);
            cto_13[i] = cto_i.getvalue(0, 2);
            cto_14[i] = cto_i.getvalue(0, 3);
            cto_15[i] = cto_i.getvalue(0, 4);
            cto_16[i] = cto_i.getvalue(0, 5);

            cto_21[i] = cto_i.getvalue(1, 0);
            cto_22[i] = cto_i.getvalue(1, 1);
            cto_23[i] = cto_i.getvalue(1, 2);
            cto_24[i] = cto_i.getvalue(1, 3);
            cto_25[i] = cto_i.getvalue(1, 4);
            cto_26[i] = cto_i.getvalue(1, 5);

            cto_31[i] = cto_i.getvalue(2, 0);
            cto_32[i] = cto_i.getvalue(2, 1);
            cto_33[i] = cto_i.getvalue(2, 2);
            cto_34[i] = cto_i.getvalue(2, 3);
            cto_35[i] = cto_i.getvalue(2, 4);
            cto_36[i] = cto_i.getvalue(2, 5);

            cto_41[i] = cto_i.getvalue(3, 0);
            cto_42[i] = cto_i.getvalue(3, 1);
            cto_43[i] = cto_i.getvalue(3, 2);
            cto_44[i] = cto_i.getvalue(3, 3);
            cto_45[i] = cto_i.getvalue(3, 4);
            cto_46[i] = cto_i.getvalue(3, 5);


            cto_51[i] = cto_i.getvalue(4, 0);
            cto_52[i] = cto_i.getvalue(4, 1);
            cto_53[i] = cto_i.getvalue(4, 2);
            cto_54[i] = cto_i.getvalue(4, 3);
            cto_55[i] = cto_i.getvalue(4, 4);
            cto_56[i] = cto_i.getvalue(4, 5);

            cto_61[i] = cto_i.getvalue(5, 0);
            cto_62[i] = cto_i.getvalue(5, 1);
            cto_63[i] = cto_i.getvalue(5, 2);
            cto_64[i] = cto_i.getvalue(5, 3);
            cto_65[i] = cto_i.getvalue(5, 4);
            cto_66[i] = cto_i.getvalue(5, 5);

        }


    }


    // writing pointers to fields

    // flag = 0  :: during nonlinear loop
    // *** dummy stress tensor for equilibrium formulation
    // *** CTO for computing Jacobian

    // flag = 1  :: after nonlinear loop
    // *** strain, stress, plastic strain, alpha, kappa
    // *** dummy stress is not needed!
    // *** CTO not needed!
    int pbox_update_params = p_update_flag[1];

    if (pbox_update_params == 1) {

        // updating mode!

        // writing new strain_tot to its field
        infields[0].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[15]);
        infields[1].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[16]);
        infields[2].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[17]);
        infields[3].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[18]);
        infields[4].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[19]);
        infields[5].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[20]);
        // writing new stress_corr to its field
        infields[6].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[21]);
        infields[7].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[22]);
        infields[8].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[23]);
        infields[9].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[24]);
        infields[10].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[25]);
        infields[11].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[26]);
        // writing new strain_p to its field
        infields[12].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[27]);
        infields[13].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[28]);
        infields[14].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[29]);
        infields[15].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[30]);
        infields[16].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[31]);
        infields[17].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[32]);
        // writing new alpha to its field
        infields[18].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[33]);
        infields[19].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[34]);
        infields[20].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[35]);
        infields[21].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[36]);
        infields[22].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[37]);
        infields[23].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[38]);
        // writing new sigma_y to its field
        infields[24].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[39]);

        // un-flagging the update
        //pbox_update_params = 0;

    } else {

        // dummy stress tensor
        infields[25].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[21]);
        infields[26].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[22]);
        infields[27].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[23]);
        infields[28].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[24]);
        infields[29].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[25]);
        infields[30].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[26]);


        // CTO
        infields[31].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[40]);
        infields[32].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[41]);
        infields[33].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[42]);
        infields[34].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[43]);
        infields[35].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[44]);
        infields[36].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[45]);

        infields[37].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[46]);
        infields[38].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[47]);
        infields[39].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[48]);
        infields[40].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[49]);
        infields[41].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[50]);
        infields[42].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[51]);

        infields[43].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[52]);
        infields[44].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[53]);
        infields[45].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[54]);
        infields[46].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[55]);
        infields[47].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[56]);
        infields[48].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[57]);

        infields[49].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[58]);
        infields[50].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[59]);
        infields[51].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[60]);
        infields[52].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[61]);
        infields[53].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[62]);
        infields[54].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[63]);

        infields[55].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[64]);
        infields[56].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[65]);
        infields[57].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[66]);
        infields[58].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[67]);
        infields[59].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[68]);
        infields[60].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[69]);

        infields[61].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[70]);
        infields[62].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[71]);
        infields[63].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[72]);
        infields[64].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[73]);
        infields[65].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[74]);
        infields[66].harmonic(1).getpointer()->setvalue(elsel, gps, meshdeform, vals[75]);

    }

    // returning strain_tot

    //return {vals[15], vals[16], vals[17],
    //        vals[18], vals[19], vals[20]};
    return {vals[15]};
}
// *****************************************************************************
// *****************************************************************************
std::vector<double> relative_stress(std::vector<double> S, std::vector<double> alpha)
{
    // it calculates the realtive stress
    // ksi = S - alpha
    std::vector<double> ksi(6);
    for (int i = 0; i < 6; i++) {
        ksi[i] = S[i] - alpha[i];
    }
    return ksi;
}
// *****************************************************************************
// *****************************************************************************
double yield_function(std::vector<double> ksi, double k)
{
    // it calculates the yield function
    // F = ||ksi|| - k = ||S - alpha|| - k
    // ksi : relative stress
    // k : kappa
    // S : deviatoric stress
    // alpha : back stress
    return (L2norm(ksi) - k);
}
// *****************************************************************************
// *****************************************************************************
std::vector<double> dev(std::vector<double> t)
{
    // it calculates the deviatoric Stress from stress
    // S = Stress - StressV
    // StressV : hydrostatic stress >> 1/3 * tr(Stress) * I

    // S : deviatoric Stress
    std::vector<double> dev_t(6);
    // p : pressure >> 1/3 * tr(Stress)
    double p = 1.0/3.0 * (t[0] + t[1] + t[2]);
    // I : eye
    std::vector<double> I{1.0, 1.0, 1.0, 0.0, 0.0, 0.0};

    for (int i = 0; i < 6; i++) {
          //
          dev_t[i] = t[i] - (p * I[i]);
    }

    return dev_t;
}
// *****************************************************************************
// *****************************************************************************
double delta_lambda(double F, double mu, double H)
{
    // it calculates delta_lambda in plastic regime
    // delta_lambda = yield_function@trial / (2*mu + H)
    // trial : the elastic predicted step
    // mu    : shear modulus
    // H     : hardening modulus


    return (F / (2.0*mu + H));
}
// *****************************************************************************
// *****************************************************************************
std::vector<double> n_hat(std::vector<double> ksi)
{
    // it calculates the unit vector (n_hat) of relative stress
    // ksi = S - alpha
    // n_hat = ksi / ||ksi||

    // L2norm of ksi
    double ksi_L2norm = L2norm(ksi);

    // declaring n_hat
    std::vector<double> nhat(6);

    // n_hat
    for (int i = 0; i < 6; i++) {
        nhat[i] = ksi[i] / ksi_L2norm;
    }

    return nhat;
}
// *****************************************************************************
// *****************************************************************************
double L2norm(std::vector<double> t)
{
    // returns the L2norm of a tensor (matrix or vector)

    // L2norm of the tensor :
    // {t_xx^2 + t_yy^2 + t_zz^2 + 2 t_yz^2 + 2 t_xz^2 + 2 t_xy^2} ^ 0.5

    double norm = 0.0;

    // the multipliers
    std::vector<double> a{1.0, 1.0, 1.0, 2.0, 2.0, 2.0};

    for (int i = 0; i < 6; i++) {
        norm = norm + ( a[i] *  t[i]*t[i] );
    }

    norm = sqrt(norm);

    return norm;
}
