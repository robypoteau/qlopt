#include <odesolver.h>

namespace thesis{

    mat OdeIntWrapper(OdeWrapper sys, const vec& x, const vec& t, double tol)
    {
        string type = "stiff";
        SUNNonlinearSolver NLS_TYPE;
        int flag; // For checking if functions have run properly
        realtype abstol = 1e-6; // real tolerance of system
        realtype reltol = 1e-6; // absolute tolerance of system

<<<<<<< HEAD
        runge_kutta4<vector_type> stepper;
        //runge_kutta_dopri5< vector_type > stepper;
        //runge_kutta_fehlberg78<vector_type> stepper;

        //bulirsch_stoer< vector_type > stepper( tol, tol, 1.0 , 1.0 );
        //bulirsch_stoer_dense_out< vector_type > stepper;

        //controlled_runge_kutta< runge_kutta_fehlberg78<vector_type> > stepper;
        //controlled_runge_kutta< runge_kutta_cash_karp54<vector_type> > stepper;
        //controlled_runge_kutta< runge_kutta_dopri5<vector_type> > stepper( 1E-6, 1E-6 );
=======
        vec ppp = sys.getParameter();

        // 3. Set problem dimensions etc. for the forward problem.
        // ---------------------------------------------------------------------
        sunindextype N_forward = x.size();
        sunindextype M = ppp.size();
        //int lt = (t(t.size()-1) - t(0))/(t(1)-t(0)) + 1;
        vec temp(N_forward*(M+1));
        mat output(N_forward*(M+1), 1);
        output.col(0) << x, vec::Zero(N_forward*M);
        // 4. Set initial conditions for the forward problem.
        // ---------------------------------------------------------------------
        N_Vector y_forward;
        N_Vector *yS0 = NULL; // Problem vector.
        y_forward = N_VNew_Serial(N_forward);
        yS0 = N_VCloneVectorArray_Serial(M, y_forward);

        for(int i = 0; i<N_forward; i++)
        {
            NV_DATA_S(y_forward)[i] = x(i);
            //cout << NV_DATA_S(y_forward)[i] << endl; //CHECKED
        }
        for (int i = 0; i < M; i++) {
            N_VConst(0.0, yS0[i]);
        }
        // 5. Create CVODES object for the forward problem.
        // ---------------------------------------------------------------------
        void *cvode_mem = NULL;
        if(type == "nonstiff"){
            NLS_TYPE = SUNNonlinSol_FixedPoint(y_forward, 0);
            cvode_mem = CVodeCreate(CV_ADAMS);
        }else if(type == "stiff"){
            NLS_TYPE = SUNNonlinSol_Newton(y_forward);
            cvode_mem = CVodeCreate(CV_BDF);
        }else{
            cerr << "You done messed up A-A-Ron! Problem must be stiff or \
                nonstiff." << endl;
            exit(0);
        }
        flag = CVodeSetNonlinearSolver(cvode_mem, NLS_TYPE);

        // ---------------------------------------------------------------------------


        auto fcvs = [](realtype t, N_Vector u, N_Vector u_dot, void *user_data){
            userData* ud = (userData*) user_data;

            auto f = [](realtype t, N_Vector u, N_Vector u_dot,
                realtype * p, OdeWrapper * sys, int N_forward){

                realtype *udata  = N_VGetArrayPointer(u);
                realtype *dudata = N_VGetArrayPointer(u_dot);
                vec temp(N_forward);

                for(int i=0; i<N_forward; i++){
                    temp(i) = (double) udata[i];
                }
                int m = sys->getParameter().size();
                vec pp(m);
                for(int i=0; i<m; i++){
                    pp(i) = (double) p[i];
                }
                temp = sys->flinear((double) t, temp, pp);

                for(int i=0; i<N_forward; i++){
                    dudata[i] = (realtype) temp(i);
                }
                return 0;
            };
            return f(t, u, u_dot, ud->p, ud->pOdeWrapper, ud->N);
        };
        //cout << ppp.transpose() << endl;
        realtype udp[M];
        for (int i = 0; i < M; i++) {
            udp[i] = (realtype) ppp(i);
        }

        userData ud = {
            .pOdeWrapper = &sys,
            .N = N_forward,
            .p = udp
        };

        flag = CVodeInit(cvode_mem, fcvs, t(0), y_forward);
        flag = CVodeSetUserData(cvode_mem, &ud);
        // 7. Specify integration tolerances for the forward problem.
        flag = CVodeSStolerances(cvode_mem, reltol, abstol);

        // 10. Create linear solver object for the forward problem.
        // ---------------------------------------------------------------------
        SUNMatrix DM = SUNDenseMatrix(N_forward, N_forward);
        SUNLinearSolver LS;
        LS = SUNLinSol_Dense(y_forward, DM);
        //LS = SUNLinSol_SPGMR(y_forward, 1, 0);
        //LS = SUNLinSol_SPFGMR(y_forward, 2, 5);
        //LS = SUNLinSol_SPBCGS(y_forward, 2, 5);
        //LS = SUNLinSol_SPTFQMR(y_forward, 0, 5);
        //LS = SUNLinSol_PCG(y_forward, 2, 5);
>>>>>>> featureCVODES


        // 12. Attach linear solver module for the forwrad problem.
        // ---------------------------------------------------------------------
        //SUNDenseMatrix J(N_forward, N_forward);
        flag = CVodeSetLinearSolver(cvode_mem, LS, DM);

        //  Initialize the cvbandpre preconditioner module
        //flag = CVBandPrecInit(cvode_mem, N_forward, 8, 8);

        //flag = CVodeSetMaxNumSteps(cvode_mem, 12001);
        //flag = CVodeSetMinStep(cvode_mem, 1e-10);

        //16 Sensitivity problem
        //----------------------------------------------------------------------
        flag = CVodeSensInit1(cvode_mem, M, CV_STAGGERED, NULL, yS0);
        flag = CVodeSensEEtolerances(cvode_mem);
        flag = CVodeSetSensParams(cvode_mem, ud.p, NULL, NULL);

        // 17. Integrate forward problem.
        // ---------------------------------------------------------------------
        //int print_steps = 120;
        auto start = high_resolution_clock::now();
        realtype tout;
        realtype step_length = t(1)-t(0);
        realtype end_time = t(t.size()-1)+step_length/2;
        realtype ts = 0;
        int ncheck = 0;
        // loop over output points, call CVode, print results, test for error
        std::cout << "Performing Forward Integration: \n\n";
        for (tout = step_length; tout <= end_time; tout += step_length) {
            flag = CVode(
                cvode_mem,
                tout,
                y_forward,
                &ts,
                CV_NORMAL);

            flag = CVodeGetSens(cvode_mem, &ts, yS0);
            for(int i = 0; i<N_forward; i++)
            {
                temp(i) = NV_DATA_S(y_forward)[i];
                //cout << NV_DATA_S(y_forward)[i] << ", "; //CHECKED
                for(int j = 0; j<M; j++)
                {
                    temp((j+1)*N_forward + i) = Ith(yS0[j],i);
                }
            }
            //cout << endl;
            output.conservativeResize(NoChange, output.cols()+1);
            output.col(output.cols()-1) = temp;
        }
        auto end = high_resolution_clock::now();
   		auto duration = duration_cast<microseconds>(end - start);
        cout << "OdeInt Time: " << duration.count()/1E6 << " s" << endl;
        //cout << output.leftCols(4) << endl;

<<<<<<< HEAD
        for( size_t i=0; i<tlen; i++ )
        {
           for(size_t j=0; j<xlen; j++ )
           {
               output(j, i) = x_vec[i][j];
           }
       }cout << output.leftCols(4)  << endl;// exit(0);
        return output;
    }
=======
        // 30. Deallocate memory.
        // ---------------------------------------------------------------------------
        N_VDestroy(y_forward);
        N_VDestroyVectorArray_Serial(yS0, M);
        CVodeFree(&cvode_mem);
        // ---------------------------------------------------------------------------
>>>>>>> featureCVODES

        // 31. Free linear solver and matrix memory for the backward problem.
        // ---------------------------------------------------------------------------
        SUNLinSolFree(LS);
        SUNNonlinSolFree(NLS_TYPE);

        return output;
    }
}
