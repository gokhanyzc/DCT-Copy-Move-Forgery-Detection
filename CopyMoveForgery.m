
clc;
clear all;

img=imread('C:\Users\Gökhan Yazýcý\Desktop\test\forged2_gj90.png');
Shift_Count=20;
Start_Time=tic;%Geçen süreyi hesaplar.
gray_img=img(:,:,1)*0.299+img(:,:,2)*0.587+img(:,:,3)*0.114; % griye çevirir.
[M N]=size(gray_img);
B=16; %Blok boyutu=B x B
MinShiftSize=B;
num_blocks=(M-B+1)*(N-B+1);%blok sayýsý
Feature_Matrix=int16(zeros(num_blocks,B*B+2));
image_mask=(zeros(256,256,3));

Q8 = ...    % standard quantization table for JPEG (8x8 matrix)
    [ 16  11  10  16  24  40  51  61; ...
    12  12  14  19  26  58  60  55; ...
    14  13  16  24  40  57  69  56; ...
    14  17  22  29  51  87  80  62; ...
    18  22  37  56  68 109 103  77; ...
    24  35  55  64  81 104 113  92; ...
    49  64  78  87 103 121 120 101; ...
    72  92  95  98 112 100 103  99 ];
Q8P=2.5*Q8;
Q8P(1,1)=2*Q8(1,1);
Q16=[Q8P                      2.5*Q8(1,8)*ones(8,8);... %quantization table for quantize 16x16 DCT blocks
    2.5*Q8(8,1)*ones(8,8)         2.5*Q8(8,8)*ones(8,8)];


%ZÝGZAG
[num_rows num_cols]=size(gray_img);

% Initialise the output vector
out=zeros(1,num_rows*num_cols);

cur_row=1;	cur_col=1;	cur_index=1;

% First element
%out(1)=in(1,1);

while cur_row<=num_rows & cur_col<=num_cols
	if cur_row==1 & mod(cur_row+cur_col,2)==0 & cur_col~=num_cols
		out(cur_index)=gray_img(cur_row,cur_col);
		cur_col=cur_col+1;							%move right at the top
		cur_index=cur_index+1;
		
	elseif cur_row==num_rows & mod(cur_row+cur_col,2)~=0 & cur_col~=num_cols
		out(cur_index)=gray_img(cur_row,cur_col);
		cur_col=cur_col+1;							%move right at the bottom
		cur_index=cur_index+1;
		
	elseif cur_col==1 & mod(cur_row+cur_col,2)~=0 & cur_row~=num_rows
		out(cur_index)=gray_img(cur_row,cur_col);
		cur_row=cur_row+1;							%move down at the left
		cur_index=cur_index+1;
		
	elseif cur_col==num_cols & mod(cur_row+cur_col,2)==0 & cur_row~=num_rows
		out(cur_index)=gray_img(cur_row,cur_col);
		cur_row=cur_row+1;							%move down at the right
		cur_index=cur_index+1;
		
	elseif cur_col~=1 & cur_row~=num_rows & mod(cur_row+cur_col,2)~=0
		out(cur_index)=gray_img(cur_row,cur_col);
		cur_row=cur_row+1;		cur_col=cur_col-1;	%move diagonally left down
		cur_index=cur_index+1;
		
	elseif cur_row~=1 & cur_col~=num_cols & mod(cur_row+cur_col,2)==0
		out(cur_index)=gray_img(cur_row,cur_col);
		cur_row=cur_row-1;		cur_col=cur_col+1;	%move diagonally right up
		cur_index=cur_index+1;
		
	elseif cur_row==num_rows & cur_col==num_cols	%obtain the bottom right element
        out(end)=gray_img(end);							%end of the operation
		break										%terminate the operation
    end
end





ShiftVectors=[0 0 0];
blockShiftIndex=[];
lastSVindx=1;

%Blocklama & Feature Extraction
rownum=1;
for i=1:M-B+1
    for j=1:N-B+1
        Block_DCT_Quantized=round(dct2(gray_img(i:i+B-1,j:j+B-1))./Q16);%DCT of Block And Quantize it
        %Bloklarý üst üste koyma
        row=zeros(1,B*B);
        for k=1:B
            row((k-1)*B+1:k*B)=Block_DCT_Quantized(k,:);
        end
        Feature_Matrix(rownum,:)=[row i j];
        rownum=rownum+1;
    end
end
%Sýralama
Feature_Matrix=sortrows(Feature_Matrix);
%Benzer bloklarý eþleþtirme
for i=1:num_blocks-1
    if(Feature_Matrix(i,1:B*B)==Feature_Matrix(i+1,1:B*B))
        shift=Feature_Matrix(i,B*B+1:B*B+2)-Feature_Matrix(i+1,B*B+1:B*B+2);
        %filtering by location difference
        if(norm(double(shift))<MinShiftSize)
            continue;
        end

        %shift vector bulma
        indx=find(ShiftVectors(:,1)==shift(1,1) & ShiftVectors(:,2)==shift(1,2));
        if(isempty(indx)==0)%found
            ShiftVectors(indx,3)=ShiftVectors(indx,3)+1;%shift counter attýr.
        else  %not found  = new shift vector
            lastSVindx=lastSVindx+1;
            ShiftVectors=[ShiftVectors; 
                             shift 1];
            indx=lastSVindx;
        end
        blockShiftIndex=[blockShiftIndex;
                            indx   i;
                            indx   i+1];
    end
end
bitis=toc(Start_Time);
fprintf('Geçen süre: %.1f Second \n', bitis);

%shift vektör blocklarýný filtreleme
for i=1:lastSVindx%shift vektörlerinin sayýsý için
    if (ShiftVectors(i,3)>Shift_Count)%shift counta göre filtreleme

        randomColor=[1,1,1];

        for j=1:size(blockShiftIndex)%Benzer blok sayýlarý için;
            if(blockShiftIndex(j,1)==i)
                %bloklarý renklendirme
                block_index=blockShiftIndex(j,2);
                block_y=Feature_Matrix(block_index,B*B+1);
                block_x=Feature_Matrix(block_index,B*B+2);
                

                for k = 0 : B-1
                    for z=0 : B-1
                        image_mask(block_y+k,block_x+z,1)=255;
                        image_mask(block_y+k,block_x+z,2)=255;
                        image_mask(block_y+k,block_x+z,3)=255;
                    end   
                     
                end
            end
        end
 
%         title(sprintf('Shift Index=%d   Shift Vector=( %d , %d )',i,ShiftVectors(i,1),ShiftVectors(i,2)));
        imshow(image_mask);
        ginput(1);
        
        maskeleme = imread('forged2_maske.png');
        test = getFmeasure(image_mask,maskeleme)
    else
        blockShiftIndex(find(blockShiftIndex(:,1)==i),:)=[];%filtrelenmiþ blocklarý kaldýrýr.
    end
end




title('Sonuç');
for i=1:size(blockShiftIndex)%benzer blok sayýsý için
                %bloklarý renklendirme
                block_index=blockShiftIndex(i,2);
                block_y=Feature_Matrix(block_index,B*B+1);
                block_x=Feature_Matrix(block_index,B*B+2);
               
                for k = 0 : B-1
                    for z=0 : B-1
                        img(block_y+k,block_x+z,1)=255;
                        img(block_y+k,block_x+z,2)=0;
                        img(block_y+k,block_x+z,3)=0;
                    end   
                end
end
imshow(img);



    