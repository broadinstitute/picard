/** reads having a soft clip in 5' larger than 2 bases  */
function accept(rec) {
	if( rec.getReadUnmappedFlag()) return false;
	var cigar = rec.getCigar();
	if( cigar == null ) return false;
	var ce = cigar.getCigarElement(0);
	return ce.getOperator().name()=="S" && ce.length()>2;
	}

accept(record);
