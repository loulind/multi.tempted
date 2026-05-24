.PHONY: clean dirs

# cleans project
clean:
	rm -rf <filepath>

# creates folders for project
dirs:
	mkdir -p <foldername>

# Creates outputs
<target1path> <target2path>: <depend1path> <depend2path>
	<script>
