function [fopt,xopt,gopt]=Gradient_F_Newton(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de gradient Newton                                        //
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

      ind = 7;
      [F,G,H] = Oracle(x,ind);
      if (G == %inf) then
          disp ('Issue with G')
       end
        
//    - test de convergence

      if norm(G) <= tol then
         kstar = k;
         break
      end

//    - calcul de la direction de descente

      if det(H)~=0 then
          disp('H inversible')
        D = -inv(H)*G
      else
          disp('H NON inversible')
          D = -G;
      end

//    - calcul de la longueur du pas de gradient

      //alpha = alphai;
      alphai = 1;
     [alpha,ok]=Wolfe(alphai,x,D,Oracle)
     

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
