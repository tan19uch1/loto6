function	x_smp	=	g_smp(b,w,N)

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
