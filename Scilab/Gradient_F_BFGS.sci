function [fopt,xopt,gopt]=Gradient_F_BFGS(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de gradient BFGS                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres du gradient a pas fixe";
   labels = ["Nombre maximal d''iterations";...
             "Valeur du pas de gradient";...
             "Seuil de convergence sur ||G||"];
   typ = list("vec",1,"vec",1,"vec",1);
   default = ["5000";"0.0005";"0.000001"];
   [ok,iter,alphai,tol] = getvalue(titre,labels,typ,default);

    exec('Wolfe_Skel.sci');

// ----------------------------
// Initialisation des variables
// ----------------------------
   
   logG = [];
   logP = [];
   Cout = [];

   timer();

// -------------------------
// Boucle sur les iterations
// -------------------------

   x = xini;

   kstar = iter;
   for k = 1:iter

//    - valeur du critere et du gradient

      ind = 4;
      if k>1 then 
          G_old = G;
          F_old = F;
      end
      [F,G] = Oracle(x,ind);
      if (G == %inf) then
          disp ('Issue with G')
       end
        
//    - test de convergence

      if norm(G) <= tol then
         kstar = k;
         break
      end

//    - calcul de la direction de descente

      if k==1 then
        W = eye(n-md,n-md);
      else
        deltax = x - x_old;
        deltag = G - G_old;
        A1 = eye(n-md,n-md)-1/(deltag'*deltax)*(deltax*deltag');
        A2 = (eye(n-md,n-md)-deltag*deltax'/(deltag'*deltax));
        C = deltax*deltax'/(deltag'*deltax);
        W = A1*W*A2+C;
        //W = (eye(n-md,n-md)-deltax*deltag'/deltag'*deltax)*W*(eye(n-md,n-md)-deltag*deltax'/deltag'*deltax)+deltax*deltax'/deltag'*deltax;
      end
      D = -W*G;

//    - calcul de la longueur du pas de gradient

      if k>1 then 
          delta = F_old - F;
          alphai = -2*delta/(G'*D);
      else
          alphai = 1;
      end
      
      
     [alpha,ok]=Wolfe(alphai,x,D,OraclePG)
     //printf('alpha: %f\n',alpha)
     //printf('ok: %i\n',ok)
     

//    - mise a jour des variables
      x_old = x;
      
      x = x + (alpha*D);

//    - evolution du gradient, du pas et du critere

      logG = [ logG ; log10(norm(G)) ];
      logP = [ logP ; log10(alpha) ];
      Cout = [ Cout ; F ];

   end

// ---------------------------
// Resultats de l'optimisation
// ---------------------------

   fopt = F;
   xopt = x;
   gopt = G;

   tcpu = timer();

   cvge = ['Iteration         : ' string(kstar);...
           'Temps CPU         : ' string(tcpu);...
           'Critere optimal   : ' string(fopt);...
           'Norme du gradient : ' string(norm(gopt))];
   disp('Fin de la methode de gradient a pas fixe')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout);

endfunction
