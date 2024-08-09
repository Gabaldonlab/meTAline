.PHONY: download-external-sources docker-image docker-container docker singularity metaline-config metaline-run metaline-debug-sandbox metaline-recompile-from-debug-sandbox

download-external-sources:
	./scripts/download-external-sources.sh

docker-image:
	@docker image build -t metaline:latest .

docker-container:
	@docker container run -dit --name metaline metaline:latest

docker:
	@docker image build -t metaline:latest .
	@docker container -dit --name metaline metaline:latest

singularity:
	./scripts/get_singularity_def_file_from_dockerfile.sh
	sudo singularity build metaline.sif metaline-singularity.def

metaline-debug-sandbox:
	sudo singularity build --sandbox ./metaline-debug metaline.sif

metaline-recompile-from-debug-sandbox:
	sudo singularity build metaline.sif ./debug
