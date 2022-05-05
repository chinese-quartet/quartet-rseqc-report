all: clean install-report install-exp2qcdt
	@echo "Compile the quartet-rseqc-report...."
	@bin/lein uberjar
	@printf "\n\n\e[1;32mRun the command for more details: \nsource .env/bin/activate\njava -jar target/uberjar/quartet-rseqc-report-*-standalone.jar -h\e[0m"

clean:
	@echo "Clean the environment..."
	@bin/lein clean
	@rm -rf .env .lsp .clj-kondo report/dist report/quartet-rseqc-report.egg-info exp2qcdt.tar.gz resources/renv/library resources/renv/staging

make-env:
	virtualenv -p python3 .env
	cp resources/bin/exp2qcdt.sh .env/bin

install-report: make-env
	. .env/bin/activate
	cd report && python3 setup.py sdist && pip3 install dist/*.tar.gz

install-exp2qcdt:
	cp -R resources/* .env/
	tar czvf .env/exp2qcdt.tar.gz exp2qcdt
	@echo 'renv::activate(".env"); renv::restore();' > .env/Rprofile
	export R_PROFILE_USER=.env/Rprofile && Rscript -e 'install.packages(".env/exp2qcdt.tar.gz", repos = NULL, type="source")'
