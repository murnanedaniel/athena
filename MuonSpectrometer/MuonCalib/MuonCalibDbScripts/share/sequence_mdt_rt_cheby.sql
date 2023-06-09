-- PROJECT: ATLAS MDT Calibration DB
-- AUTHOR: Elisabetta.Vilucchi@lnf.infn.it, celio@roma3.infn.it
-- DATE: May 08
-- VERSION v1r1

--USAGE:
--sqlplus <USER_NAME>/<USER_PASSWORD> @sequence_mdt_rt_cheby.sql <USER_NAME>

--DROP SEQUENCE &1.."MDT_RT_CHEBY_SEQUENCE";
CREATE SEQUENCE &1.."MDT_RT_CHEBY_SEQUENCE" NOCYCLE NOORDER NOCACHE NOMAXVALUE MINVALUE 1 INCREMENT BY 1 START WITH 1;
