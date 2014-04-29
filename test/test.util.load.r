# Test extraction of samples by metadata
"test.extract.samples.by.metadata" <- function(input="./../code/test/test.extract.samples.by.metadata.input.txt"){
    extract.samples.by.metadata(input, 'GENDER:M;COUNTRY:!USA')
    extract.samples.by.metadata(input, 'GENDER:*;COUNTRY:!USA')
    extract.samples.by.metadata(input, 'GENDER:*;COUNTRY:*,!USA')
    extract.samples.by.metadata(input, 'GENDER:*;COUNTRY:*,!USA;DISEASE:CD')
    extract.samples.by.metadata(input, 'GENDER:*;COUNTRY:*,!USA;DISEASE:UC')
    extract.samples.by.metadata(input, 'GENDER:MALE;COUNTRY:*,!USA;DISEASE:UC')
    extract.samples.by.metadata(input, 'GENDER:M;COUNTRY:*,!USA;DISEASE:UC')
}