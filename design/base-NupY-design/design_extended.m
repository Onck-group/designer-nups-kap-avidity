%makeNup3X.m
clear all


%%%%%%%%% Turn making a protein ON (=1) or OFF (=0)
makeSeq=1;
A_or_B='B';


%Data from AF:

%vector of AA's in order, together with the AA densities of both (in
%percentage) in the SPACER contents in both domains, resp.
if A_or_B=='A'
    aaVector =         ['R'   'K'     'D'    'E'    'N' 'Q' 'S' 'T' 'H' 'C' 'Y' 'G' 'A' 'P' 'M' 'V' 'W' 'I' 'L' 'F'];
    collAADensVector = [0.4,   2.7,   0.0,   0.1,   18.4, 10.5, 16.3, 13.8,  0.0,  0.1,  0.2, 13.9,  8.0,  5.1,  1.5,  1.4,  0.1,  1.3,  3.2,  3.1]/100;
    extAADensVector =  [0.0370, 0.1073, 0.0767, 0.0657, 0.0937, 0.0347, 0.1450, 0.0577, 0.0120, 0.0000, 0.0123, 0.0263, 0.0697, 0.0527, 0.0210, 0.0297, 0.0033, 0.0503, 0.0767, 0.0297];									
else
    aaVector =         ['R'       'K'       'D'       'E'       'N'       'Q'       'S'       'T'       'H'       'C'       'Y'       'G'       'A'       'P'       'M'       'V'       'W'       'I'       'L'       'F'];
    collAADensVector = [0.003151, 0.025835, 0.000630, 0.001260, 0.196597, 0.119723, 0.168872, 0.137366, 0.000000, 0.000630, 0.001890, 0.124134, 0.077505, 0.049779, 0.013863, 0.015753, 0.000630, 0.012602, 0.035287, 0.014493];
    extAADensVector  = [0.036395, 0.107452, 0.076256, 0.065858, 0.093588, 0.034662, 0.143847, 0.057192, 0.012132, 0.000000, 0.013865, 0.025997, 0.069324, 0.051993, 0.020797, 0.029463, 0.003466, 0.050260, 0.077990, 0.029463];
end

seqNup49=['MFGLNKASSTPAGGLFGQASGASTGNANTGFSFGGTQTGQNTGPSTGGLFGAKPAGSTGGLGASFGQQQQQSQTNAFGGSATTGGGLFGNKPNNTANTGGGLFGANSNSNSGSLFGSNNAQTSRGLFGNNNTNNINNSSSGMNNASAGLFGSKPAGGTSLFGNTSTSSAPAQNQGMFGAKPAGTSLFGNNAGNTTTGGGLFGSKPTGATSLFGSSNNNNNNNNSNNIMSASGGLFGNQQQQLQQQPQMQCALQNLSQLPITPMTRISELPPQIRQEIEQLDQYIQKQVQISHHLKADTIDHDELIDSIPRDVAYLLKSESATSQYLKQDLKKISSFKSLIDEDLLDTQTFSVLLQQLLTPGSKISSNDLDKFFQKKIHLYEKKLEDYCRILSDIETAVNGIDTDLFGAPNNPNSTAITADLGSSEAENLLQLKTGLAAIVSTVIEEFTLFMDIAERIAVLHQKTKTLASLSI'];
seqNup49C=seqNup49(1:251);
seqNup49E=[''];

seqNup57=['MFGFSGSNNGFGNKPAGSTGFSFGQNNNNTNTQPSASGFGFGGSQPNSGTATTGGFGANQATNTFGSNQQSSTGGGLFGNKPALGSLGSSSTTASGTTATGTGLFGQQTAQPQQSTIGGGLFGNKPTTTTGGLFGNSAQNNSTTSGGLFGNKVGSTGSLMGGNSTQNTSNMNAGGLFGAKPQNTTATTGGLFGSKPQGSTTNGGLFGSGTQNNNTLGGGGLFGQSQQPQTNTAPGLGNTVSTQPSFAWSKPSTGSNLQQQQQQQIQVPLQQTQAIAQQQQLSNYPQQIQEQVLKCKESWDPNTTKTKLRAFVYNKVNETEAILYTKPGHVLQEEWDQAMEKKPSPQTIPIQIYGFEGLNQRNQVQTENVAQARIILNHILEKSTQLQQKHELDTASRILKAQSRNVEIEKRILKLGTQLATLKNRGLPLGIAEEKMWSQFQTLLQRSEDPAGLGKTNELWARLAILKERAKNISSQLDSKLMVFNDDTKNQDSMSKGTGEESNDRINKIVEILTNQQRGITYLNEVLEKDAAIVKKYKNKT'];
seqNup57C=seqNup57(1:255);
seqNup57E=[''];

seqNup100=['MFGNNRPMFGGSNLSFGSNTSSFGGQQSQQPNSLFGNSNNNNNSTSNNAQSGFGGFTSAAGSNSNSLFGNNNTQNNGAFGQSMGATQNSPFGSLNSSNASNGNTFGGSSSMGSFGGNTNNAFNNNSNSTNSPFGFNKPNTGGTLFGSQNNNSAGTSSLFGGQSTSTTGTFGNTGSSFGTGLNGNGSNIFGAGNNSQSNTTGSLFGNQQSSAFGTNNQQGSLFGQQSQNTNNAFGNQNQLGGSSFGSKPVGSGSLFGQSNNTLGNTTNNRNGLFGQMNSSNQGSSNSGLFGQNSMNSSTQGVFGQNNNQMQINGNNNNSLFGKANTFSNSASGGLFGQNNQQQGSGLFGQNSQTSGSSGLFGQNNQKQPNTFTQSNTGIGLFGQNNNQQQQSTGLFGAKPAGTTGSLFGGNSSTQPNSLFGTTNVPTSNTQSQQGNSLFGATKLTNMPFGGNPTANQSGSGNSLFGTKPASTTGSLFGNNTASTTVPSTNGLFGNNANNSTSTTNTGLFGAKPDSQSKPALGGGLFGNSNSNSSTIGQNKPVFGGTTQNTGLFGATGTNSSAVGSTGKLFGQNNNTLNVGTQNVPPVNNTTQNALLGTTAVPSLQQAPVTNEQLFSKISIPNSITNPVKATTSKVNADMKRNSSLTSAYRLAPKPLFAPSSNGDAKFQKWGKTLERSDRGSSTSNSITDPESSYLNSNDLLFDPDRRYLKHLVIKNNKNLNVINHNDDEASKVKLVTFTTESASKDDQASSSIAASKLTEKAHSPQTDLKDDHDESTPDPQSKSPNGSTSIPMIENEKISSKVPGLLSNDVTFFKNNYYISPSIETLGNKSLIELRKINNLVIGHRNYGKVEFLEPVDLLNTPLDTLCGDLVTFGPKSCSIYENCSIKPEKGEGINVRCRVTLYSCFPIDKETRKPIKNITHPLLKRSIAKLKENPVYKFESYDPVTGTYSYTIDHPVLT'];
seqNup100C=seqNup100(1:610);
seqNup100E=seqNup100(611:800);

seqNup116=['MFGVSRGAFPSATTQPFGSTGSTFGGQQQQQQPVANTSAFGLSQQTNTTQAPAFGNFGNQTSNSPFGMSGSTTANGTPFGQSQLTNNNASGSIFGGMGNNTALSAGSASVVPNSTAGTSIKPFTTFEEKDPTTGVINVFQSITCMPEYRNFSFEELRFQDYQAGRKFGTSQNGTGTTFNNPQGTTNTGFGIMGNNNSTTSATTGGLFGQKPATGMFGTGTGSGGGFGSGATNSTGLFGSSTNLSGNSAFGANKPATSGGLFGNTTNNPTNGTNNTGLFGQQNSNTNGGLFGQQQNSFGANNVSNGGAFGQVNRGAFPQQQTQQGSGGIFGQSNANANGGAFGQQQGTGALFGAKPASGGLFGQSAGSKAFGMNTNPTGTTGGLFGQTNQQQSGGGLFGQQQNSNAGGLFGQNNQSQNQSGLFGQQNSSNAFGQPQQQGGLFGSKPAGGLFGQQQGASTFASGNAQNNSIFGQNNQQQQSTGGLFGQQNNQSQSQPGGLFGQTNQNNNQPFGQNGLQQPQQNNSLFGAKPTGFGNTSLFSNSTTNQSNGISGNNLQQQSGGLFQNKQQPASGGLFGSKPSNTVGGGLFGNNQVANQNNPASTSGGLFGSKPATGSLFGGTNSTAPNASSGGIFGSNNASNTAATTNSTGLFGNKPVGAGASTSAGGLFGNNNNSSLNNSNGSTGLFGSNNTSQSTNAGGLFQNNTSTNTSGGGLFSQPSQSMAQSQNALQQQQQQQRLQIQNNNPYGTNELFSKATVTNTVSYPIQPSATKIKADERKKASLTNAYKMIPKTLFTAKLKTNNSVMDKAQIKVDPKLSISIDKKNNQIAISNQQEENLDESILKASELLFNPDKRSFKNLINNRKMLIASEEKNNGSQNNDMNFKSKSEEQETILGKPKMDEKETANGGERMVLSSKNDGEDSATKHHSRNMDEENKENVADLQKQEYSEDDKKAVFADVAEKDASFINENYYISPSLDTLSSYSLLQLRKVPHLVVGHKSYGKIEFLEPVDLAGIPLTSLGGVIITFEPKTCIIYANLPNRPKRGEGINVRARITCFNCYPVDKSTRKPIKDPNHQLVKRHIERLKKNPNSKFESYDADSGTYVFIVNHAAEQT'];
seqNup116C=seqNup116(172:764);
seqNup116E=seqNup116(765:960);

seqNup145N=['MFNKSVNSGFTFGNQNTSTPTSTPAQPSSSLQFPQKSTGLFGNVNVNANTSTPSPSGGLFNANSNANSISQQPANNSLFGNKPAQPSGGLFGATNNTTSKSAGSLFGNNNATANSTGSTGLFSGSNNIASSTQNGGLFGNSNNNNITSTTQNGGLFGKPTTTPAGAGGLFGNSSSTNSTTGLFGSNNTQSSTGIFGQKPGASTTGGLFGNNGASFPRSGETTGTMSTNPYGINISNVPMAVADMPRSITSSLSDVNGKSDAEPKPIENRRTYSFSSSVSGNAPLPLASQSSLVSRLSTRLKATQKSTSPNEIFSPSYSKPWLNGAGSAPLVDDFFSSKMTSLAPNENSIFPQNGFNFLSSQRADLTELRKLKIDSNRSAAKKLKLLSGTPAITKKHMQDEQDSSENEPIANADSVTNIDRKENRDNNLDNTYLNGKEQSNNLNKQDGENTLQHEKSSSFGYWCSPSPEQLERLSLKQLAAVSNFVIGRRGYGCITFQHDVDLTAFTKSFREELFGKIVIFRSSKTVEVYPDEATKPMIGHGLNVPAIITLENVYPVDKKTKKPMKDTTKFAEFQVFDRKLRSMREMNYISYNPFGGTWTFKVNHFSIWGLVNEEDAEIDEDDLSKQEDGGEQPLRKVRTLAQSKPSDKEVILKTDGTFGTLSGKDDSIVEEKAYEPDLSDADFEGIEASPKLDVSKDWVEQLILAGSSLRSVFATSKEFDGPCQNEIDLLFSECNDEIDNAKLIMKERRFTASYTFAKFSTGSMLLTKDIVGKSGVSIKRLPTELQRKFLFDDVYLDKEIEKVTIEARKSNPYPQISESSLLFKDALDYMEKTSSDYNLWKLSSILFDPVSYPYKTDNDQVKMALLKKERHCRLTSWIVSQIGPEIEEKIRNSSNEIEQIFLYLLLNDVVRASKLAIESKNGHLSVLISYLGSNDPRIRDLAELQLQKWSTGGCSIDKNISKIYKLLSGSPFEGLFSLKELESEFSWLCLLNLTLCYGQIDEYSLESLVQSHLDKFSLPYDDPIGVIFQLYAANENTEKLYKEVRQRTNALDVQFCWYLIQTLRFNGTRVFSKETSDEATFAFAAQLEFAQLHGHSLFVSCFLNDDKAAEDTIKRLVMREITLLRASTNDHILNRLKIPSQLIFNAQALKDRYEGNYLSEVQNLLLGSSYDLAEMAIVTSLGPRLLLSNNPVQNNELKTLREILNEFPDSERDKWSVSINVFEVYLKLVLDNVETQETIDSLISGMKIFYDQYKHCREVAACCNVMSQEIVSKILEKNNPSIGDSKAKLLELPLGQPEKAYLRGEFAQDLMKCTYKI'];
seqNup145NC=seqNup145N(1:242);
seqNup145NE=seqNup145N(243:433);





seqColl=  seqNup100C;           %=seqNup100(1:610);
seqExt=   seqNup100E;          %seqNup100(611:800);

%remove GLFGs from coll. seq
[seqCollGLFG,seqCollnoGLFG]= regexp(seqColl,'GLFG','match','split');
originalNumberGLFG=length(seqCollGLFG);
seqCollnoGLFG=strjoin(seqCollnoGLFG,'');
%repeat to remove left-over FG's
[seqCollFGnoGLFG,seqCollnoFGnoGLFG]=regexp(seqCollnoGLFG,'FG','match','split');
seqCollnoFGnoGLFG=strjoin(seqCollnoFGnoGLFG,'');
originalNumberFG=length(seqCollFGnoGLFG);





for indexProtein = 1:50
    if makeSeq==1
        %set domain lengths of Nup100 unfolded domain here
        collDomainLength        =   length(seqColl);
        extDomainLength         =   length(seqExt);

        %10+FG+10+GLFG repeat unit = 26 length
        numPatternRepeats = floor(collDomainLength/111);
        numFGRepeats = floor(collDomainLength/111)*4+2;
        numGLFGRepeats = floor(collDomainLength/111)*3+1;
        numRepeats=numFGRepeats+numGLFGRepeats;


        %set domain lengths of Nup100 unfolded domains - FG/GLFG occurences here
        collSpacerDomainLength  =   collDomainLength-numFGRepeats*2-numGLFGRepeats*4;
        extSpacerDomainLength   =   extDomainLength;

        
        occurencesCollSpacer    =   round(collSpacerDomainLength * collAADensVector) ;
        occurencesExtSpacer   =   round(extSpacerDomainLength * extAADensVector) ;
        for i = 1:length(aaVector)
            fprintf('aa is: %s,\t occ is: %f, unrounded: %f\n',aaVector(i),occurencesExtSpacer(i),extSpacerDomainLength*extAADensVector(i))
            
        end
  
        





        % Below pattern contains 4:3 FG, GLFG in a 90-sequence repeat pattern:
        % 10_FG_10_GLFG_10_FG_10_GLFG_10_FG_10_GLFG_10_FG
        %fgSeq=['          FG          GLFG          FG          GLFG          FG          GLFG          FG          FG          GLFG          FG          GLFG          FG          GLFG          FG          FG          GLFG          FG          GLFG          FG          GLFG          FG          FG          GLFG          FG          GLFG          FG          GLFG          FG          FG          GLFG          FG          GLFG          FG          GLFG          FG          FG          GLFG          FG          GLFG          FG          GLFG          FG          FG          GLFG          FG          GLFG          FG      '];
        fgSeq=['             FG             GLFG             FG             GLFG             FG             GLFG             FG             FG             GLFG             FG             GLFG             FG             GLFG             FG             FG             GLFG             FG             GLFG             FG             GLFG             FG             FG             GLFG             FG             GLFG             FG             GLFG             FG             FG             GLFG             FG             GLFG             FG             GLFG             FG             FG             GLFG             FG        '];

        collSpacerSeq='';
        extSpacerSeq='';

        for i=1:length(aaVector)
            currentAA=aaVector(i);
            for j=1:occurencesCollSpacer(i)
                collSpacerSeq=strjoin({collSpacerSeq,aaVector(i)},'');
            end
            for k=1:occurencesExtSpacer(i)
                extSpacerSeq=strjoin({extSpacerSeq,aaVector(i)},'');
            end
        end

        %fill up with Glycine if the length of a domain is not met d.t. rounding
        while length(extSpacerSeq) < extSpacerDomainLength
                extSpacerSeq=strjoin({extSpacerSeq,'S'},'');
                fprintf('ping!')
        end

        while length(collSpacerSeq) < collSpacerDomainLength
                collSpacerSeq=strjoin({collSpacerSEq,'S'},'');
                
        end


        %create the extended domain and spacer seqs for the collapsed
        %domain. Keep repeating these loop until no accidental FG or GF
        %motifs are present. Note that we do this separately for both
        %domains (since the chances of both domains having nothing is
        %small, more efficient now).
        
        no_accidental_FG_GF = false ;
        while no_accidental_FG_GF == false
            
            %re-seed random sequence, initialize empty string array
            numRandColl=randsample(collSpacerDomainLength,collSpacerDomainLength);
            collSpacerSeqRand='';
            for i=1:length(numRandColl)
                collSpacerSeqRand=strjoin({collSpacerSeqRand,collSpacerSeq(numRandColl(i))},'');
            end
            
            %check whether FG and GLFG exist:
            contains_FG_motif=regexp(collSpacerSeqRand,'FG');
            contains_GF_motif=regexp(collSpacerSeqRand,'GF');
            if isempty(contains_FG_motif) && isempty(contains_GF_motif)
                no_accidental_FG_GF = true;
            end
        end
        
        no_accidental_FG_GF = false ;
        while no_accidental_FG_GF == false 
            numRandExt=randsample(extSpacerDomainLength,extSpacerDomainLength);
            extSpacerSeqRand='';
            for i=1:length(numRandExt)
                extSpacerSeqRand=strjoin({extSpacerSeqRand,extSpacerSeq(numRandExt(i))},'');
            end
            
            contains_FG_motif=regexp(extSpacerSeqRand,'FG');
            contains_GF_motif=regexp(extSpacerSeqRand,'GF');
            if isempty(contains_FG_motif) && isempty(contains_GF_motif)
                no_accidental_FG_GF = true;
                extSeq = extSpacerSeqRand ; % this is the one we pass on
                
            end
            
        end
        


        %make the collapsed domain
        fgCounter=0;
        glfgCounter=0;
        patternCounter=0;

        for i=1:length(collSpacerSeqRand)
            fgSeq(1,i+fgCounter*2+glfgCounter*4)=collSpacerSeqRand(numRandColl(i));


            if mod(patternCounter,2)==0
                if mod(i,13)==0 && mod(i,26)~=0
                    fgCounter=fgCounter+1;
                end
                if mod(i,26)==0
                    glfgCounter=glfgCounter+1;
                end
            end

            if mod(patternCounter,2)==1
                if mod(i,13)==0 && mod(i,26)~=0
                    glfgCounter=glfgCounter+1;
                end
                if mod(i,26)==0
                    fgCounter=fgCounter+1;
                end
            end

            if mod(i,91)==0
                patternCounter=patternCounter+1;
            end
        end



        finalSeq=strjoin({fgSeq,extSeq},'');
        
        %Write the full protein sequence
       % protName=['NupY_' num2str(indexProtein) '.txt'];
       % fid=fopen(['./output_sequences/' protName],'w');
       % fprintf(fid,'%%NupY_%d\n',indexProtein);
       % fprintf(fid,'%s',finalSeq);
       % fclose(fid);
        
        %Write the collapsed domain
       % protNameC=['NupYC_' num2str(indexProtein) '.txt'];
       % fid=fopen(['./output_sequences/' protNameC],'w');
       % fprintf(fid,'%%NupYC_%d\n',indexProtein);
       % fprintf(fid,'%s',fgSeq);
       % fclose(fid);
        
        %Write the extended domain
        protNameE=['NupYE_' num2str(indexProtein) '.txt'];
        fid=fopen(['./only_ext_domains/' protNameE],'w');
        %fprintf(fid,'%%NupYE_%d\n',indexProtein);
        fprintf(fid,'%s',extSeq);
        fclose(fid);
    end
end


    