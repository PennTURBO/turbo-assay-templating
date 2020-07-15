```sql
SELECT
	PK_LAB_RESULT_ITEM_ID ,
	LAB_RESULT_ITEM_CODE ,
	LAB_RESULT_ITEM_DESCRIPTION ,
	LOINC
FROM
	mdm.R_LAB_RESULT_ITEM rlri
WHERE
	lower(LAB_RESULT_ITEM_DESCRIPTION) LIKE '%sars%cov%2%'
	AND LOINC IS NOT NULL
```



| PK_LAB_RESULT_ITEM_ID | LAB_RESULT_ITEM_CODE   | LAB_RESULT_ITEM_DESCRIPTION | LOINC   | PDS mentions |
| --------------------- | ---------------------- | --------------------------- | ------- | -----------: |
| 20669003              | SARS-COV-2 IGG         | SARS-COV-2 IGG              | 94563-4 |          621 |
| 20529020              | SARS COV 2 RNA, RT PCR | SARS COV 2 RNA, RT PCR      | 94309-2 |            2 |

