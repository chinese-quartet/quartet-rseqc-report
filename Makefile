all: clean make-env install-report install-exp2qcdt
	@echo "Compile the quartet-rseqc-report...."
	@bin/lein uberjar
	@printf "\n\n\e[1;32mRun the command for more details: \nsource .env/bin/activate\njava -jar target/uberjar/quartet-rseqc-report-*-standalone.jar -h\e[0m"

clean:
	@echo "Clean the environment..."
	@bin/lein clean
	@rm -rf .env .lsp .clj-kondo report/dist report/quartet-rseqc-report.egg-info exp2qcdt.tar.gz resources/renv/library resources/renv/staging

make-env:
	pip3 install virtualenv && virtualenv -p python3 .env
	cp resources/bin/exp2qcdt.sh .env/bin

install-report:
	cd report && .env/bin/python3 setup.py sdist && .env/bin/pip3 install dist/*.tar.gz

install-exp2qcdt:
	@Rscript -e 'install.packages("renv", repos="http://cran.us.r-project.org")'
	cp -R resources/* .env/
	@echo 'renv::activate(".env"); renv::restore();' > .env/Rprofile
	export R_PROFILE_USER=.env/Rprofile && Rscript -e 'renv::install("./exp2qcdt");'
