context("test extract_sites")

site = c(2, 5, 9, 10, 15, 20)

test_that("test extract_sites", {

	expect_that(extract_sites(0, 1, site, index = FALSE), equals(integer(0)))
	expect_that(extract_sites(0, 2, site, index = FALSE), equals(2))
	expect_that(extract_sites(1, 8, site, index = FALSE), equals(c(2, 5)))
	expect_that(extract_sites(2, 5, site, index = FALSE), equals(c(2, 5)))
	expect_that(extract_sites(3, 4, site, index = FALSE), equals(integer(0)))
	expect_that(extract_sites(1, 21, site, index = FALSE), equals(c(2, 5, 9, 10, 15, 20)))
	expect_that(extract_sites(9, 16, site, index = FALSE), equals(c(9, 10, 15)))
	expect_that(extract_sites(15, 20, site, index = FALSE), equals(c(15, 20)))
	expect_that(extract_sites(20, 25, site, index = FALSE), equals(c(20)))
	expect_that(extract_sites(25, 26, site, index = FALSE), equals(integer(0)))
})

if(Sys.getenv("IS_PBS") != "") {

site = sort(sample(10000000, 1000000))
pos = do.call("rbind", lapply(1:1000, function(i) sort(sample(max(site), 2))))

system.time(for(i in 1:1000) {
	which(site >= pos[i, 1] & site <= pos[i, 2])
})

site2 = as.double(site)
pos2 = as.double(pos)
dim(pos2) = dim(pos)
system.time(for(i in 1:1000) {
	extract_sites(pos2[i, 1], pos2[i, 2], site2)
})

}
