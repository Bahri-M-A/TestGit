function test_loubna()

% Algorithme servant au traitement des sections transversales
%
% Cet algorithme est utilis� pour l?analyse des coupes transversales 
% dans le but de d�terminer le diam�tre �quivalent des sections 
% et de caract�riser de la texture de celles-ci.
% Ces op�rations n�cessitent la binarisation des images.
% Une analyse des niveaux de gris est ensuite effectu�e de mani�re 
% � �valuer les profils d?humidit�.   

% =========================================================================
% Effacement des donn�es pr�c�dentes 
clear all 
% Cr�ation et ouvertures des diff�rents fichiers de stockage 
fid = fopen ('reswet.xls', 'a' ) ; 
fida = fopen ('resdenswet.xls', 'a' ) ; 
fidb = fopen ('resdensvwet.xls', 'a' ) ; 
fidc = fopen ('resintsolprowet.xls', 'a' ) ; 
fidd=fopen('reswetcrevglobal.xls','a'); 
fide=fopen('reswetcrevarea.xls','a'); 
fidf=fopen('reswetcreveccen.xls','a'); 
% Impression d'une ligne de commentaire dans le premier fichier de stockage 
fprintf(fid,'Aire (mm�) \t Facteur forme \t excentricit� \t Diam�tre �quivalent (mm) \tAire crevasse (mm�) \t Facteur de forme \t Excentricit� crevasses \t Intensit� moyenne \t Ecart-type \t ratio crevasse \t Intensit� moyenne solide \n') ; 
 
%Boucle permettant d'analyse l'ensemble des sections 
 for jj=107:10:347,
     
        % Chargement du fichier '.dat' obtenu apr�s reconstruction des donn�es tomographiques 

        s=['C:\Angel\these\Retinne\tomo\13-02-01\wet\A_0' num2str(jj) '.dat']; 
        s1=load(s); 
        taille=size(s1); 
        m=sqrt(taille(1,1)); 
        ss=reshape(s1,m,m); 
        % Passage du fichier '.dat' � un fichier image '.tif' 
        imwrite(ss,['C:\Angel\these\Retinne\tomo\13-02-01\A_0' num2str(jj) .tif'], 'tif','compression', 'none'); 
        delete(s);   clear ss  
        file=strcat('A_0',num2str(jj),'.tif'); 
        i=imread(file); 
        % D�finition de la r�solution 
        scale=41.0e-3; 


        % D�termination du seuil de binarisation par la m�thode d'Otsu 
        level = graythresh(i); 
        % Binarisation 
        bin = im2bw(i,level); 
        % Erosion pour �liminer les artefacts issus du seuillage 
        er1=bwmorph(bin,'erode',1); 
        % Reconstruction 
        fin = imreconstruct(er1,bin); 
        % Obtention de l'image pleine 

        % S�rie de 10 fermetures 
        clo=bwmorph(fin,'dilate',10); 
        ero=bwmorph(clo,'erode',10); 
        % Remplissage de l'image 
        fill1= imfill(ero,'holes'); 
        % Ouverture pour adoucir le contour 
        fill=bwmorph(fill1,'open',1); 


        % Obtention de l'image des crevasses 

        % Ou exclusif entre l'image finale et l'image apr�s binarisation 
        crev1=xor(fin,fill); 
        % Labellisation de l'image pour d�terminer les diff�rents objets qui la composent 
        [crev2,NUM] = BWLABEL(crev1,8); 
        % S�lection des crevasses dont l'aire est sup�rieure � 10 pixels 
        STATS = REGIONPROPS(crev2,'Area'); 
        idx=find([STATS.Area]>10); 
        crev3=ismember(crev2,idx); 
        % Labellisation de l'image finale utilis�e pour caract�riser les crevasses 
        [crev4,NUM] = BWLABEL(crev3,8); 

        % Mesure de l'aire et de l'excentricit� de chaque crevasse 

         STATS = REGIONPROPS(crev4,'Area','Eccentricity'); 

        All1=[[STATS.Area]];  
         All2=[[STATS.Eccentricity] ]; 

        % Impression des r�sultats     
        fprintf(fidd,'image %d\n',jj); 
        fprintf(fidd,'\r'); 
        fprintf(fidd,'%f\t',All1) ; 
        fprintf(fidd,'\r'); 
        fprintf(fidd,'%f\t',All2) ; 
        fprintf(fidd,'\r'); 
        fprintf(fide,'%f\t',All1) ; 
        fprintf(fide,'\r');
        fprintf(fidf,'%f\t',All2) ; 
        fprintf(fidf,'\r'); 

        % Effacement des donn�es 
        clear All1;clear All2;clear All3; 

        % Mesures effectu�es de mani�re globale sur l'image pleine et sur l'image des crevasses 

        % Mesure de l'aire de la section 
        surfp1=bwarea(fill);surfp=(scale^2)*surfp1; 
        % Mesure du p�rim�tre de la section 
        perp1=bwperim(fill);perp=scale*bwarea(perp1); 
        % Calcul du facteur de forme de la section 
        shapp=(perp^2)/(4.*pi*surfp);  
        % Mesure de l'aire des crevasses 
        surfc1=bwarea(crev3);surfc=(scale^2)*surfc1; 
        % Mesure du p�rim�tre des crevasses 
        perc1=bwperim(crev3);perc=scale*bwarea(perc1); 
        % Calcul du facteur de forme des crevasses 
        shapc=(perc^2)/(4.*pi*surfc); 
        % Calcul du diam�tre �quivalent de la section 
        dequi=sqrt(4*surfp/pi); 
        % Calcul du ration crevasse 
        ratiocrev=surfc/surfp; 
        eccplein=regionprops(fill,'Eccentricity'); 
        % Calcul de l'excentricit� de la section remplie 
        eccpl=[eccplein.Eccentricity]; 
        % Calcul de l'excentricit� des crevasses 
        ecccrev=regionprops(crev3,'eccentricity'); 
        ecccr=[ecccrev.Eccentricity]; 
        % Calcul du niveau de gris moyen de l'image et de son �cart-type 
        grey=immultiply(i,fill); 
        moyenne=sum(sum(grey))/surfp1; 
        m=size(grey); 
        M=(moyenne*ones(m)); 
        N=minus(double(grey),M); 
        variance=(sum(sum(N.*N))-sum(sum(1-double(fill)))*(moyenne)^2)/(surfp1-1); 
        ecarttype=sqrt(variance); 
        % Calcul de l'intensit� ramen�e � la fraction du mat�riau humide 
        intsol=moyenne/(1-ratiocrev); 
        % Impression des r�sultats dans le fichier de stockage 
        fprintf(fid,'%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n',surfp,shapp,eccpl,dequi,surfc,shapc,ecccr,moyenne,ecarttype,ratiocrev,intsol) ; 
        clear i;clear bin;clear er1; clear clo;clear fill1;clear fin;clear crev1; 


        % Analyse en niveau de gris 

        % Calcul du nombre de pseudo anneaux successifs 
        no=not(fill); 
        distance=bwdist(no); 
        lim=max(max(distance))/sqrt(2); 

        % Construction du premier pseudo anneau 
        pas(1)=5*scale;  
        filla=fill; 
        fillb=bwmorph(filla,'erode',5); 
        aro=xor(filla,fillb); 

        % Mesure de l'aire de l'image originale masqu�e par le pseudo anneau occup�e par des 
        %crevasses 
        resa=and(aro,crev3); 
        res=bwarea(resa);  
        % Mesure de l'aire du pseudo anneau 
        norm=bwarea(aro);  
        % Calcul du ratio de crevasse dans le pseudo anneau de l'image originale 
        dens(1)=res/norm;  
        % Calcul du niveau de gris moyen de l'image originale masqu�e par le pseudo anneau  
        multi=immultiply(grey,aro); 
        volume=sum(sum(multi)); 
        densv(1)=volume/norm; 
        % Calcul du niveau de gris rapport� � la fraction de mat�riau humide 
        intsolpro(1)=densv(1)/(1-dens(1)); 

        % Construction des anneaux suivants et mesures des diff�rentes grandeurs 

        k=2; 
        for m=5:5:lim 

            pas(k)=(m+5)*scale;  
            filla=bwmorph(fill,'erode',m); 
            fillb=bwmorph(filla,'erode',5); 
            aro=xor(filla,fillb); 

            resa=and(aro,crev3); 
            res=bwarea(resa);  
            norm=bwarea(aro);  
            dens(k)=res/norm;  

            multi=immultiply(grey,aro); 
            volume=sum(sum(multi)); 
            densv(k)=volume/norm; 
            intsolpro(k)=densv(k)/(1-dens(k));   

            k=k+1; 
        end 
        fin=k-1; 

        % Stokage des donn�es dans les fichiers    

        resdens=[dens]; 
        resdensv=[densv]; 
        resintsolpro=[intsolpro]; 
        fprintf(fida,'%f\t',resdens) ; 
        fprintf(fida,'\r'); 
        fprintf(fidb,'%f\t',resdensv) ; 
        fprintf(fidb,'\r'); 
        fprintf(fidc,'%f\t',resintsolpro) ; 
        fprintf(fidc,'\r'); 
        clear resdens;clear resdensv;clear resintsolpro;clear dens;clear densv;clear 
        intsolpro,clear pas 

 end % end boucle for
 
% Effacement des donn�es 
clear aro;clear crev;clear fill; clear filla,clear fillb; 
clear resa;clear norm;clear multi;clear volume; 
% Fermeture des fichiers de stockage 
fclose(fid); 
fclose(fida); 
fclose(fidb); 
fclose(fidc); 
fclose(fidd); 
fclose(fide); 
fclose(fidf);

end %end fuction