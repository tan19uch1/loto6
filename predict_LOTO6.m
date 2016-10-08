function	predict_LOTO6()
	
	%-------------------------------------------
	%
	%	Predicting ROTO6
	%				using Boltzman Machine
	%
	%		Author:Shotaro Taniguchi
	%
	%-------------------------------------------
	
	%parameter
	N_smp	=	500;
	N_lrn	=	100;
	num		=	43;
	eta_b	=	0.0025;
	eta_w	=	0.0005;
	
	%initialize
	b	=	rand(1,num);
	w	=	rand(num,num);
	w	=	w	-	diag(diag(w));

	%data load
	if	exist('./loto6.mat')
		fprintf(1,'LOADING ... \n');
		load './loto6.mat'
	else
		fprintf(1,'MAKING TEST DATA ... \n');
		loto6_csv		=	csvread('./loto6.csv');
		loto6_number	=	loto6_csv(2:end,4:9);
		x	=	[];
		for	row	=	loto6_number'
			x_row	=	zeros(1,43);
			x_row(row)	=	1;
			x	=	[x;x_row];
		end
		testdata	=	x;
		save('loto6.mat','x','testdata')
	end

	%learning
	fprintf(1,'LEARNING ... \n');
	histdata	=	repmat([1:43],size(x,1),1);
	histdata	=	reshape(histdata',1,[]);
	x_hist_test	=	reshape(x',1,[]);
	for	n	=	1:1:N_lrn
		fprintf(1,'*--iter = %d\n',n);
		fflush(stdout);
		%renew parameter (bi)
		x_smp	=	gibbs_sampling(b,w,N_smp);
		L_b		=	sum(x)	-	sum(x_smp).*(size(x,1)./size(x_smp,1));
		b		=	b	+	eta_b.*L_b;

		%renew parameter (wij)
		for	i	=	1:1:num
			for	j	=	1:1:num
				L_wij	=	sum(x(:,i).*x(:,j))	-	sum(x_smp(:,i).*x_smp(:,j)).*(size(x,1)./size(x_smp,1));
				w(i,j)	=	w(i,j)	+	eta_w.*L_wij;
		 	end
		end
		w	=	w	-	diag(diag(w));
		save('./result.mat')
		%see histgram
		x_hist	=	reshape(x_smp',1,[]);
		f	=	figure(1);
		title('HISTOGRAM')
			f1	=	subplot(2,1,1);
			hist(histdata(find(x_hist)),43);
			xlim([1,43])
			title(f1,'SAMPLING DATA by BM')
			xlabel('NUMBER (1-43)')
			ylabel('FREQUENCY')
			f2	=	subplot(2,1,2);
			hist(histdata(find(x_hist_test)),43);
			xlim([1,43])
			title(f1,'TEST DATA')
			xlabel('NUMBER (1-43)')
			ylabel('FREQUENCY')
	end

	
end


%------------------
%	Gibbs Sampling
%------------------
function	x_smp	=	gibbs_sampling(b,w,N)

	%parameter
	burn	=	50;
	wait	=	5;
	cnt		=	0;

	%initialize
	x		=	randn(size(b))>1;
	x_smp	=	[];
	
	%monte carlo sampling
	while	size(x_smp,1)	<	N+burn
		psb_vec	=	[];
		%sampling
		for	i_smp	=	1:1:length(b)
%			i_smp
			psb	=	p_i(i_smp,b,w,x);
			pu	=	rand;
			if	psb >=	pu
				x(i_smp)	=	1;
			else
				x(i_smp)	=	0;
			end
		end
		%judge
		if	rem(cnt,wait)==0;
			x_smp	=	[x_smp;x];
		end
		cnt	=	cnt + 1;
	end

	%burn-in
	fprintf(1,'\n');
	x_smp	=	x_smp(end-N+1:end,:);
	
end

function	psb_i	=	p_i(i_smp,b,w,x)

	x_tmp			=	x;
	x_tmp(i_smp)	=	0;
	lambda		=	b(i_smp)	+	sum(x_tmp*w./2);
	psb_i		=	exp(lambda)./(1+exp(lambda));

end
