numRuns = 5;
for pBase = [ 500, 1000, 3000,  5000, 7000,  9000, 10000, 20000]
	for s=[0.05,0.10,  0.15]
		for oBase = [2, 5, 10]
			sBase = ceil(s*pBase);
			n = ceil(oBase*sBase*log(pBase));
			for eps = [0.05, 0.1, 0.5, 1,1.75]
    				samples = genData(n,pBase,sBase,eBase,eps,numRuns);
			    	filename = sprintf('/data/Samples_p%d_s%d_e%g_o%g_Ceps%g',pBase,sBase,eBase,oBase,eps);    
			    	filenameM = sprintf('Samples_p%d_s%d_e%g_o%g_Ceps%g.mat',pBase,sBase,eBase,oBase,eps);
			    	mkdir (filename)
				filename
			   	pBase
				sBase
				oBase
				cd (filename)
				save(filenameM,'samples');
			    	theta = samples.theta;
				X = samples.X;
				y = samples.y;
				save('theta','theta','-ascii');
				save('X','X','-ascii');
			       	save('y','y','-ascii');
		        	cd '~/sparse-recovery/code/IHT_Code'
			end
		end
	end
end
