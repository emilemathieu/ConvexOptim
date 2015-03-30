function [F,G,ind]=OraclePG(qc,ind)
    //qc : vecteur des d ́ebits r ́eduits,
    //F : valeur du crit`ere  ́evalu ́e au point qC ,
    //G : vecteur des d ́eriv ́ees du crit`ere par rapport `a qC,
    //ind : indicateur du type de calcul `a effectuer : 
        //= 2 : calcul de F uniquement,
        //= 3 : calcul de G uniquement,
        //=4: calculdeF etG.
        
        //F = 0;
        //G = 0;
        
        if ind==2 then
            F = 1/3*sum((q0+B*qc).*(r.*(q0+B*qc).*abs(q0+B*qc)))+sum(pr.*(Ar*(q0+B*qc)));
        elseif ind==3 then
            G =B'*(r.*abs(q0+B*qc).*(q0+B*qc))+B'*Ar'*pr; 
        elseif ind==4 then
            a = q0+B*qc;
            F = 1/3*sum(a.*(r.*a.*abs(a)))+sum(pr.*(Ar*a));
            //F = 1/3*sum((q0+B*qc).*(r.*(q0+B*qc).*abs(q0+B*qc)))+sum(pr.*(Ar*(q0+B*qc)));
            G =B'*(r.*abs(q0+B*qc).*(q0+B*qc))+B'*Ar'*pr;
        end

endfunction
