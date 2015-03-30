function [alphan,ok]=Wolfe(alpha,x,D,Oracle)


//////////////////////////////////////////////////////////////
//                                                          //
//   RECHERCHE LINEAIRE SUIVANT LES CONDITIONS DE WOLFE     //
//                                                          //
//                                                          //
//  Arguments en entree                                     //
//  -------------------                                     //
//    alpha  : valeur initiale du pas                       //
//    x      : valeur initiale des variables                //
//    D      : direction de descente                        //
//    Oracle : nom de la fonction Oracle                    //
//                                                          //
//  Arguments en sortie                                     //
//  -------------------                                     //
//    alphan : valeur du pas apres recherche lineaire       //
//    ok     : indicateur de reussite de la recherche       //
//             = 1 : conditions de Wolfe verifiees          //
//             = 2 : indistinguabilite des iteres           //
//                                                          //
//                                                          //
//    omega1 : coefficient pour la 1-ere condition de Wolfe //
//    omega2 : coefficient pour la 2-eme condition de Wolfe //
//                                                          //
//////////////////////////////////////////////////////////////


// -------------------------------------
// Coefficients de la recherche lineaire
// -------------------------------------

   omega1 = 0.1;
   omega2 = 0.9;

   alphamin = 0.0;
   alphamax = %inf;

   ok = 0;
   dltx = 0.00000001;
   
   itermax = 1000000;
   iter = 0;

// ---------------------------------
// Algorithme de Fletcher-Lemarechal
// ---------------------------------

   // Appel de l'oracle au point initial
   
   ind = 4;
   [F,G] = Oracle(x,ind);

   // Initialisation de l'algorithme

   alphan = alpha;
   xn     = x;

   // Boucle de calcul du pas
   //
   // xn represente le point pour la valeur courante du pas,
   // xp represente le point pour la valeur precedente du pas.

//disp('-----------------------------------------------------------------')

   while ok == 0
      
      xp = xn;
      xn = x + (alphan*D);
      [Fp,Gp] = Oracle(xp,ind);
      [Fn,Gn] = Oracle(xn,ind);

      // Calcul des conditions de Wolfe
      //printf('iter is %i\n',iter)
      if Fn>(F+omega1*alphan*G'*D) then
          //disp('first Wolf condition NOT verified')
          alphamax = alphan;
          alphan = 0.5*(alphamin+alphamax);
      elseif Gn'*D<omega2*G'*D then
          //disp('second Wolf condition NOT verified')
          alphamin = alphan;
          if alphamax == %inf then
              alphan = 2*alphamin;
          else
              alphan = 0.5*(alphamin+alphamax);
          end
      else
          //disp('BOTH Wolf condition verified')
          ok = 1;
      end

      // Test d'indistinguabilite

      if norm(xn-xp) < dltx then
        ok = 2;
      end
      
      // Test nombre d'itérations
      iter = iter + 1;
      if iter >= itermax then
        ok = 3;
      end
      //printf('alphan: %f\n',alphan)
      //printf('alphamin: %f\n',alphamin)
      //printf('alphamax: %f\n',alphamax)
   end


endfunction
