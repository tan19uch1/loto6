function	predict_ROTO_SIX()
	
	%-------------------------------------------
	%
	%	predicting ROTO6
	%				using Boltzman Machine
	%
	%		author:Shotaro Taniguchi
	%
	%-------------------------------------------
	
	%parameter
	N_smp	=	50;
	N_lrn	=	100;
	num		=	43;
	eta_b	=	0.05;
	eta_w	=	0.01;
	
	%initialize
	b	=	rand(1,num);
	w	=	rand(num,num);
	w	=	w	-	diag(diag(w));

	%data load
	x	=	zeros(100,num);
	x(:,18)	=	rand(100,1)>0.5;
	x(:,19)	=	rand(100,1)>0.2;
	x(:,20)	=	1;
	x(:,21)	=	rand(100,1)>0.1;
	x(:,22)	=	rand(100,1)>0.4;
	x(:,23)	=	rand(100,1)>0.8;

	%learning
	fprintf(1,'LEARNING ... \n');
	for	n	=	1:1:N_lrn
	x	=	zeros(100,num);
	x(:,18)	=	rand(100,1)>0.5;
	x(:,19)	=	rand(100,1)>0.2;
	x(:,20)	=	1;
	x(:,21)	=	rand(100,1)>0.1;
	x(:,22)	=	rand(100,1)>0.4;
	x(:,23)	=	rand(100,1)>0.8;
		fprintf(1,'*--iter = %d\n',n);
		fflush(stdout);
		%renew parameter (bi)
		x_smp	=	gibbs_sampling(b,w,N_smp);
[sum(x_smp).*(size(x,1)./size(x_smp,1));sum(x)]
norm(sum(x_smp)-sum(x))
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
	end

	
end


%------------------
%	Gibbs Sampling
%------------------
function	x_smp	=	gibbs_sampling(b,w,N)

	%parameter
	burn	=	5;

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
		if	sum(x)	==	6 || true;
			x_smp	=	[x_smp;x];
			fprintf(1,'.');
			fflush(stdout);
		end
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
